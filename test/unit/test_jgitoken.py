#! /usr/bin/env python3
"""Offline tests for finding a usable JGI session token (`mycotools.lib.jgitoken`).

Nothing here touches the real clipboard, environment or terminal: each test
starts from "nothing available" and adds exactly the source it is about.
"""

import sys
import json
import time
import base64
import subprocess

import pytest

from mycotools.lib import jgitoken
from mycotools.lib.jgitoken import (
    MIN_SECONDS_LEFT,
    describe_left,
    is_jgi_token,
    is_usable,
    jwt_claims,
    read_clipboard,
    resolve_jgi_token,
    seconds_left,
    token_from_environment,
)

JGI_ISS = "https://auth.jgi.doe.gov/auth/realms/jgi-web"


def make_jwt(left=3600, iss=JGI_ISS, **extra):
    """A token shaped like the Data Portal's: header.payload.signature, with the
    payload's `exp` `left` seconds from now."""

    def enc(obj):
        raw = json.dumps(obj).encode()
        return base64.urlsafe_b64encode(raw).rstrip(b"=").decode()

    claims = {"iss": iss, "exp": time.time() + left, **extra}
    return f"{enc({'alg': 'RS256'})}.{enc(claims)}.signature"


class FakeStdin:
    def __init__(self, tty):
        self._tty = tty

    def isatty(self):
        return self._tty


@pytest.fixture(autouse=True)
def nothing_available(monkeypatch):
    """No clipboard, no environment, no terminal unless a test adds one."""
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: None)
    monkeypatch.delenv(jgitoken.ENV_TOKEN, raising=False)
    monkeypatch.delenv(jgitoken.ENV_TOKEN_FILE, raising=False)
    monkeypatch.setattr(sys, "stdin", FakeStdin(False))


# --------------------------------------------------------------------------- #
# reading a token's own expiry
# --------------------------------------------------------------------------- #
def test_expiry_is_read_without_the_network():
    left = seconds_left(make_jwt(left=7200))
    assert 7190 < left <= 7200


def test_expired_token_has_negative_time_left():
    assert seconds_left(make_jwt(left=-60)) < 0
    assert describe_left(make_jwt(left=-60)) == "expired"


@pytest.mark.parametrize("junk", ["", "abc", "a.b.c", "a.b", None, 12345])
def test_non_jwt_has_no_readable_expiry(junk):
    assert jwt_claims(junk) is None
    assert seconds_left(junk) is None


def test_usable_needs_time_left_but_trusts_an_opaque_token():
    assert is_usable(make_jwt(left=7200))
    assert not is_usable(make_jwt(left=-1))
    # a token about to run out would not survive one genome's download
    assert not is_usable(make_jwt(left=MIN_SECONDS_LEFT - 30))
    assert is_usable("opaque-token-with-no-readable-expiry")
    assert not is_usable("")
    assert not is_usable(None)


def test_describe_left_formats_hours_and_minutes():
    assert describe_left(make_jwt(left=11 * 3600 + 52 * 60 + 30)) == "11h52m"
    assert describe_left(make_jwt(left=25 * 60 + 30)) == "25m"
    assert describe_left("opaque") == "unknown time"


def test_only_a_jgi_issued_jwt_counts_as_a_jgi_token():
    assert is_jgi_token(make_jwt())
    assert not is_jgi_token(make_jwt(iss="https://accounts.example.com/"))
    assert not is_jgi_token("just some text that was on the clipboard")


# --------------------------------------------------------------------------- #
# the clipboard
# --------------------------------------------------------------------------- #
def _fake_run(stdout=b"", returncode=0, raises=None):
    def run(command, **kwargs):
        if raises:
            raise raises
        return subprocess.CompletedProcess(command, returncode, stdout=stdout)

    return run


def test_clipboard_text_is_returned(monkeypatch):
    # the autouse fixture stubs the reader out; test the real one here
    monkeypatch.undo()
    monkeypatch.setattr(subprocess, "run", _fake_run(b"Bearer abc"))
    assert read_clipboard() == "Bearer abc"


@pytest.mark.parametrize(
    "run",
    [
        _fake_run(raises=FileNotFoundError()),  # no pbpaste/xclip installed
        _fake_run(raises=subprocess.TimeoutExpired("x", 3)),  # hung (no display)
        _fake_run(b"", returncode=0),  # empty clipboard
        _fake_run(b"junk", returncode=1),  # tool reports failure
    ],
)
def test_unreadable_clipboard_is_just_nothing(monkeypatch, run):
    monkeypatch.undo()
    monkeypatch.setattr(subprocess, "run", run)
    assert read_clipboard() is None


def test_huge_clipboard_is_truncated(monkeypatch):
    monkeypatch.undo()
    monkeypatch.setattr(subprocess, "run", _fake_run(b"x" * 10**6))
    assert len(read_clipboard()) == jgitoken._MAX_CLIPBOARD_BYTES


# --------------------------------------------------------------------------- #
# environment
# --------------------------------------------------------------------------- #
def test_token_from_environment_variable():
    assert token_from_environment({jgitoken.ENV_TOKEN: "Bearer abc"}) == "abc"


def test_token_from_file_is_reread_each_call(tmp_path):
    path = tmp_path / "tok"
    env = {jgitoken.ENV_TOKEN_FILE: str(path)}
    path.write_text("first\n")
    assert token_from_environment(env) == "first"
    path.write_text("second\n")  # a long unattended run is handed a fresh token
    assert token_from_environment(env) == "second"


def test_unreadable_token_file_is_nothing(tmp_path):
    env = {jgitoken.ENV_TOKEN_FILE: str(tmp_path / "missing")}
    assert token_from_environment(env) is None
    assert token_from_environment({}) is None


# --------------------------------------------------------------------------- #
# resolve_jgi_token: the order things are tried in
# --------------------------------------------------------------------------- #
def test_nothing_anywhere_gives_none():
    assert resolve_jgi_token() is None


def test_explicit_token_wins_over_the_clipboard(monkeypatch):
    explicit, clip = make_jwt(), make_jwt()
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: clip)
    assert resolve_jgi_token(explicit=explicit) == explicit


def test_clipboard_token_is_used_with_no_flags_at_all(monkeypatch):
    """The whole point: copy the token, run the command."""
    clip = make_jwt()
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: "Bearer " + clip)
    assert resolve_jgi_token() == clip


def test_unrelated_clipboard_text_is_never_used(monkeypatch):
    for text in ("my shopping list", make_jwt(iss="https://other.example.com/")):
        monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3, t=text: t)
        assert resolve_jgi_token() is None


def test_expired_clipboard_token_is_not_used(monkeypatch):
    monkeypatch.setattr(
        jgitoken, "read_clipboard", lambda timeout=3: make_jwt(left=-3600)
    )
    assert resolve_jgi_token() is None


def test_expired_explicit_token_falls_through_to_the_clipboard(monkeypatch):
    fresh = make_jwt()
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: fresh)
    assert resolve_jgi_token(explicit=make_jwt(left=-3600)) == fresh


def test_environment_beats_the_clipboard(monkeypatch):
    env_token, clip = make_jwt(), make_jwt()
    monkeypatch.setenv(jgitoken.ENV_TOKEN, env_token)
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: clip)
    assert resolve_jgi_token() == env_token


def test_stored_token_is_used_only_while_unexpired():
    stored = make_jwt()
    assert resolve_jgi_token(stored=stored) == stored
    # months-old stored token (updates are infrequent): ignored, not returned
    assert resolve_jgi_token(stored=make_jwt(left=-86400 * 60)) is None


def test_a_rejected_token_is_not_offered_again(monkeypatch):
    rejected = make_jwt()  # unexpired by its own clock, but JGI said no
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: rejected)
    assert resolve_jgi_token(exclude=rejected) is None
    other = make_jwt()
    assert resolve_jgi_token(explicit=other, exclude=rejected) == other


# --------------------------------------------------------------------------- #
# the prompt
# --------------------------------------------------------------------------- #
def _answers(monkeypatch, *typed):
    """Script what the user types at the hidden prompt; returns the call log."""
    asked = []
    it = iter(typed)

    def getpass(prompt=""):
        asked.append(prompt)
        return next(it)

    monkeypatch.setattr(jgitoken.getpass, "getpass", getpass)
    return asked


def test_no_prompt_without_a_terminal(monkeypatch):
    asked = _answers(monkeypatch, "should never be read")
    assert resolve_jgi_token() is None
    assert asked == []  # a cron job / CHTC node must not hang waiting for input


def test_prompt_enter_reads_the_clipboard(monkeypatch):
    """Expired mid-run: the user copies a fresh token and presses Enter."""
    monkeypatch.setattr(sys, "stdin", FakeStdin(True))
    fresh = make_jwt()
    # read once by the silent clipboard check (nothing yet), once when Enter is
    # pressed after the user has copied a fresh token
    clipboards = iter([None, fresh])
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: next(clipboards))
    asked = _answers(monkeypatch, "")
    assert resolve_jgi_token(reason="expired") == fresh
    assert len(asked) == 1


def test_prompt_accepts_a_pasted_token(monkeypatch):
    monkeypatch.setattr(sys, "stdin", FakeStdin(True))
    pasted = make_jwt()
    _answers(monkeypatch, "Bearer " + pasted)
    assert resolve_jgi_token() == pasted


def test_prompt_rejects_an_expired_paste_then_takes_a_good_one(monkeypatch):
    monkeypatch.setattr(sys, "stdin", FakeStdin(True))
    good = make_jwt()
    _answers(monkeypatch, make_jwt(left=-100), good)
    assert resolve_jgi_token() == good


def test_prompt_gives_up_after_three_tries(monkeypatch):
    monkeypatch.setattr(sys, "stdin", FakeStdin(True))
    asked = _answers(monkeypatch, "", "", "")
    assert resolve_jgi_token() is None
    assert len(asked) == 3


def test_prompt_ctrl_d_gives_up_cleanly(monkeypatch):
    monkeypatch.setattr(sys, "stdin", FakeStdin(True))

    def eof(prompt=""):
        raise EOFError

    monkeypatch.setattr(jgitoken.getpass, "getpass", eof)
    assert resolve_jgi_token() is None


def test_the_token_is_never_logged(monkeypatch, caplog):
    secret = make_jwt()
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: secret)
    with caplog.at_level("DEBUG"):
        assert resolve_jgi_token() == secret
    assert secret not in caplog.text
    assert secret.split(".")[1] not in caplog.text


# --------------------------------------------------------------------------- #
# `mtdb update`: the same token handling, with no JGI prompt of its own
# --------------------------------------------------------------------------- #
@pytest.fixture
def update_mod(monkeypatch):
    from mycotools.mtdb import update

    calls = []
    monkeypatch.setattr(
        update, "login_check", lambda **kw: calls.append(kw) or (None, None)
    )
    update._login_calls = calls
    return update


def test_update_asks_for_no_jgi_token_when_mycocosm_is_off(update_mod, monkeypatch):
    monkeypatch.setattr(
        update_mod, "resolve_jgi_token", lambda **kw: pytest.fail("asked for a token")
    )
    assert update_mod.gather_login(False) == (None, None)


def test_update_takes_the_token_from_the_clipboard(update_mod, monkeypatch):
    clip = make_jwt()
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: "Bearer " + clip)
    assert update_mod.gather_login(True) == (None, clip)


def test_update_never_lets_login_check_prompt_for_the_jgi_token(update_mod, monkeypatch):
    """login_check's own prompt would ask for a token that may be 2 months stale
    before the clipboard was even looked at."""
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: make_jwt())
    update_mod.gather_login(True)
    assert update_mod._login_calls == [{"jgi": False}]


def test_update_ignores_a_stored_token_months_out_of_date(update_mod, monkeypatch):
    fresh = make_jwt()
    monkeypatch.setattr(
        update_mod, "login_check", lambda **kw: (None, make_jwt(left=-86400 * 60))
    )
    monkeypatch.setattr(jgitoken, "read_clipboard", lambda timeout=3: fresh)
    assert update_mod.gather_login(True) == (None, fresh)


def test_update_exits_cleanly_when_mycocosm_is_wanted_and_no_token_exists(update_mod):
    with pytest.raises(SystemExit) as exc:
        update_mod.gather_login(True)
    assert exc.value.code == 100
