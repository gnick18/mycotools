#! /usr/bin/env python3
"""
Find a usable JGI Data Portal session token with as little effort from the
user as possible.

JGI only issues these tokens to a person logged in through ORCID with two-factor
authentication, so nothing can mint one unattended. They last 12 hours. The
Data Portal's "Copy My Session Token" puts one on the clipboard, so the
workflow this module supports is: copy the token, run the command, done. When a
long run outlives its token the same thing happens again - copy a fresh one and
press Enter.

A token is looked for, in order, in

    1. the value passed explicitly (``--token`` / ``--jgi_token``),
    2. the ``JGI_SESSION_TOKEN`` environment variable, or the file named by
       ``JGI_SESSION_TOKEN_FILE`` (re-read on every call, so a long unattended
       run can be given a fresh token by overwriting the file),
    3. the clipboard,
    4. the credential store (``mtdb manage --store``),
    5. an interactive prompt (hidden input, or Enter to read the clipboard).

A candidate is only used while it still has time left. The expiry is read from
the token itself (it is a JWT) without any network call, and nothing is ever
logged except how long the chosen token remains valid. Tokens are not verified
here; the server's answer is the authority (``jgi.jgi_api_login``).
"""

import os
import sys
import json
import time
import base64
import getpass
import logging
import subprocess

from mycotools.lib.dbtools import clean_jgi_token

logger = logging.getLogger(__name__)

ENV_TOKEN = "JGI_SESSION_TOKEN"
ENV_TOKEN_FILE = "JGI_SESSION_TOKEN_FILE"

# a token with less than this left is treated as expired: it would not outlast
# even one genome's download
MIN_SECONDS_LEFT = 600
# the clipboard can hold anything; a token is ~1.5 kB
_MAX_CLIPBOARD_BYTES = 16384
_JGI_HOST = "jgi.doe.gov"
_PROMPT_ATTEMPTS = 3


def jwt_claims(token):
    """Decode a JWT's payload WITHOUT verifying it. Returns a dict, or None if
    `token` is not a decodable JWT."""
    try:
        parts = str(token).split(".")
        if len(parts) != 3:
            return None
        payload = parts[1] + "=" * (-len(parts[1]) % 4)
        claims = json.loads(base64.urlsafe_b64decode(payload))
    except (ValueError, TypeError):
        return None
    return claims if isinstance(claims, dict) else None


def seconds_left(token, now=None):
    """Seconds until the token expires (negative once it has), or None when the
    token carries no readable expiry."""
    claims = jwt_claims(token)
    if not claims or not isinstance(claims.get("exp"), (int, float)):
        return None
    return claims["exp"] - (time.time() if now is None else now)


def is_jgi_token(token):
    """True if `token` is a JWT issued by JGI's login service. Used to decide
    whether text found on the clipboard is a token at all."""
    claims = jwt_claims(token)
    return bool(claims) and _JGI_HOST in str(claims.get("iss", ""))


def is_usable(token, min_left=MIN_SECONDS_LEFT, now=None):
    """True unless the token is known to have run out. A token with no readable
    expiry is given the benefit of the doubt; the server is the authority."""
    if not token:
        return False
    left = seconds_left(token, now)
    return left is None or left >= min_left


def describe_left(token, now=None):
    """Human-readable time remaining, e.g. '11h52m', for log lines."""
    left = seconds_left(token, now)
    if left is None:
        return "unknown time"
    if left <= 0:
        return "expired"
    hours, rest = divmod(int(left), 3600)
    minutes = rest // 60
    return f"{hours}h{minutes:02d}m" if hours else f"{minutes}m"


def read_clipboard(timeout=3):
    """Return the clipboard's text, or None if it cannot be read. Never raises
    and never prompts: every failure (no clipboard tool, no display, an SSH
    session) just means 'nothing found'."""
    if sys.platform == "darwin":
        commands = [["pbpaste"]]
    elif sys.platform.startswith("win"):
        commands = [["powershell", "-NoProfile", "-Command", "Get-Clipboard"]]
    else:
        commands = [
            ["wl-paste", "--no-newline"],
            ["xclip", "-selection", "clipboard", "-o"],
            ["xsel", "--clipboard", "--output"],
        ]
    for command in commands:
        try:
            out = subprocess.run(
                command,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                timeout=timeout,
                check=False,
            )
        except (OSError, subprocess.SubprocessError):
            continue
        if out.returncode == 0 and out.stdout:
            return out.stdout[:_MAX_CLIPBOARD_BYTES].decode("utf-8", "ignore")
    return None


def token_from_environment(environ=None):
    """The token named by JGI_SESSION_TOKEN, or held in the file named by
    JGI_SESSION_TOKEN_FILE. Returns None if neither is set or readable."""
    environ = os.environ if environ is None else environ
    if environ.get(ENV_TOKEN):
        return clean_jgi_token(environ[ENV_TOKEN])
    path = environ.get(ENV_TOKEN_FILE)
    if path:
        try:
            with open(os.path.expanduser(path), "r") as raw:
                return clean_jgi_token(raw.read())
        except OSError:
            logger.warning(f"Cannot read {ENV_TOKEN_FILE}={path}")
    return None


def _prompt(reason):
    """Ask the user for a token. Empty input reads the clipboard, so the whole
    interaction is: copy the token, press Enter."""
    print(flush=True)
    if reason:
        print(reason, flush=True)
    print(
        "Log in at https://data.jgi.doe.gov, open the avatar menu (top right) and "
        "click 'Copy My Session Token'.",
        flush=True,
    )
    for _ in range(_PROMPT_ATTEMPTS):
        try:
            typed = getpass.getpass(
                prompt="Press Enter to use the token on your clipboard, or paste it here: "
            )
        except (EOFError, KeyboardInterrupt):
            return None
        if typed.strip():
            candidate = clean_jgi_token(typed)
        else:
            # the clipboard can hold anything, so only a JGI token counts
            clip = clean_jgi_token(read_clipboard())
            candidate = clip if clip and is_jgi_token(clip) else None
        if candidate and is_usable(candidate):
            return candidate
        if candidate:
            logger.warning(
                f"That token is not usable ({describe_left(candidate)}); copy a "
                "fresh one from the Data Portal"
            )
        else:
            logger.warning("No JGI session token found on the clipboard")
    return None


def resolve_jgi_token(
    explicit=None, stored=None, prompt=True, reason=None, exclude=None
):
    """Return a usable JGI session token, or None if none could be found.

    `explicit` is a token the caller already holds and `stored` one from the
    credential store; both are used only while unexpired. `exclude` is a token
    that was just rejected by the server and must not be offered again. `prompt`
    allows the interactive fallback, which is only ever attempted on a terminal.
    `reason` heads that prompt (e.g. "The JGI session token has expired.")."""

    def acceptable(token):
        return bool(token) and token != exclude and is_usable(token)

    explicit = clean_jgi_token(explicit)
    if acceptable(explicit):
        return explicit
    if explicit and explicit != exclude:
        logger.warning(
            f"The supplied JGI session token is not usable "
            f"({describe_left(explicit)}); looking for another"
        )

    token = token_from_environment()
    if acceptable(token):
        logger.info(f"Using the JGI session token from the environment "
                    f"(valid for {describe_left(token)})")
        return token

    clip = clean_jgi_token(read_clipboard())
    if clip and is_jgi_token(clip) and acceptable(clip):
        logger.info(f"Using the JGI session token on your clipboard "
                    f"(valid for {describe_left(clip)})")
        return clip

    stored = clean_jgi_token(stored)
    if acceptable(stored):
        return stored
    if stored and stored != exclude:
        logger.info(
            f"The stored JGI session token is not usable ({describe_left(stored)}); "
            "ignoring it"
        )

    if prompt and sys.stdin.isatty():
        return _prompt(reason)
    return None
