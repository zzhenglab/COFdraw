"""Optional provider transport is mocked; no real API credentials or paid calls."""
from copy import deepcopy
from email.message import Message
from io import BytesIO
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import HTTPHandler, Request, build_opener

import pytest

from cof_vision import cofdraw_llm as llm


KEY = "unit-test-secret-4815162342"
URL = "https://provider.example/v1"
LIBRARY = {
    "sources": [{"id": "source-one", "recipe_id": "source-one", "label": "Source One", "source_registered": True,
                 "qualifications": ["Finite provisional source drawing"], "private_project": "do not send"}],
    "monomers": [{"id": "aldehyde", "label": "Trialdehyde", "full_name": "Named trialdehyde", "smiles": "DO_NOT_SEND_SMILES"},
                 {"id": "amine", "label": "Diamine", "compatible_roles": [{"linkage": "imine", "role": "primary_arylamine", "reactive_sites": 2}]}],
    "linkages": [{"id": "imine", "label": "Imine", "supported": True},
                 {"id": "unsupported_link", "label": "Reference only", "supported": False}],
    "templates": [{"id": "hcb_3_2", "label": "Hexagonal", "instances": 12, "allow_pairs": [[3, 2], [2, 3]]}],
    "styles": [{"id": "source", "label": "Source style"}, {"id": "standard", "label": "Standard"}],
    "palettes": [{"id": "classic", "label": "Classic"}],
    "saved_projects": {"secret_project": "NEVER_SEND_PRIVATE_DATA"},
}
PROPOSAL = {"selection": {"source_recipe_id": "source-one", "style_id": "source"},
            "explanation": "The named source matches the description.", "unresolved": []}


class Response(BytesIO):
    def __init__(self, data, *, endpoint=URL + "/chat/completions", headers=None, code=200):
        super().__init__(data)
        self.endpoint, self.code = endpoint, code
        self.headers = Message()
        for key, value in (headers or {}).items():
            self.headers[key] = value

    def geturl(self):
        return self.endpoint

    def getcode(self):
        return self.code

    def info(self):
        return self.headers


def wire(monkeypatch, content=None, *, raw=None, error=None, **response_kwargs):
    if raw is None:
        content = json.dumps(PROPOSAL if content is None else content) if not isinstance(content, str) else content
        raw = json.dumps({"choices": [{"finish_reason": "stop", "message": {"role": "assistant", "content": content}}]}).encode()
    calls = []
    class Opener:
        def open(self, request, timeout):
            calls.append((request, timeout))
            if error:
                raise error
            return Response(raw, **response_kwargs)
    def factory(*handlers):
        assert any(isinstance(h, llm._NoRedirect) for h in handlers)
        assert any(isinstance(h, llm.ProxyHandler) and h.proxies == {} for h in handlers)
        return Opener()
    monkeypatch.setattr(llm, "build_opener", factory)
    return calls


def select(**kwargs):
    return llm.select_with_llm("Use Source One with its source style", LIBRARY,
        base_url=kwargs.pop("base_url", URL), model=kwargs.pop("model", "chosen-provider-model"),
        api_key=kwargs.pop("api_key", KEY), **kwargs)


def test_success_has_one_call_key_only_in_authorization_and_minimal_catalog(monkeypatch, capsys):
    original = deepcopy(LIBRARY)
    calls = wire(monkeypatch)
    assert select() == PROPOSAL
    assert LIBRARY == original
    assert len(calls) == 1
    request, timeout = calls[0]
    assert request.full_url == URL + "/chat/completions"
    assert request.get_header("Authorization") == "Bearer " + KEY
    assert timeout == llm.SOCKET_TIMEOUT
    body = json.loads(request.data)
    assert set(body) == {"model", "messages", "stream"}
    assert body["model"] == "chosen-provider-model" and body["stream"] is False
    assert KEY.encode() not in request.data
    assert b"DO_NOT_SEND_SMILES" not in request.data
    assert b"NEVER_SEND_PRIVATE_DATA" not in request.data
    assert b"do not send" not in request.data
    assert b"unsupported_link" not in request.data
    assert capsys.readouterr() == ("", "")
    prompt = json.loads(body["messages"][1]["content"])
    role = next(r for r in prompt["catalog"]["monomers"] if r["id"] == "amine")["compatible_roles"][0]
    assert role == {"linkage": "imine", "role": "primary_arylamine", "reactive_sites": 2}


def test_valid_design_ids_and_unresolved_proposal_are_returned_without_execution(monkeypatch):
    proposal = {"selection": {"monomer_a_id": "aldehyde", "monomer_b_id": "amine", "linkage": "imine",
                             "template_id": "hcb_3_2", "style_id": "standard", "palette": "classic"},
                "explanation": "Catalog candidates only.", "unresolved": ["Confirm intended source identity."]}
    wire(monkeypatch, proposal)
    assert select() == proposal


@pytest.mark.parametrize("content", ["not json", "```json\n{}\n```", '{"selection":{},"selection":{},"explanation":"","unresolved":[]}',
    {**PROPOSAL, "recipe": {}}, {**PROPOSAL, "selection": {"smiles": "C"}},
    {**PROPOSAL, "selection": {"source_recipe_id": "invented"}},
    {**PROPOSAL, "selection": {"linkage": "unsupported_link"}},
    {**PROPOSAL, "selection": {"style_id": None}},
    {**PROPOSAL, "selection": {"source_recipe_id": "source-one", "monomer_a_id": "aldehyde"}},
    {**PROPOSAL, "explanation": "x" * 601}, {**PROPOSAL, "unresolved": ["x"] * 9},
    {**PROPOSAL, "unresolved": ["\x00bad"]}])
def test_invalid_or_invented_output_is_rejected_without_echo(monkeypatch, content):
    wire(monkeypatch, content)
    with pytest.raises(llm.LLMSelectionError) as caught:
        select()
    assert KEY not in str(caught.value)
    assert URL not in str(caught.value)
    assert "invented" not in str(caught.value)


@pytest.mark.parametrize("envelope", [
    {"choices": []}, {"choices": [{"message": {"content": "{}"}}] * 2},
    {"choices": [{"finish_reason": "length", "message": {"content": json.dumps(PROPOSAL)}}]},
    {"choices": [{"message": {"content": json.dumps(PROPOSAL), "tool_calls": [{"name": "run"}]}}]},
    {"choices": [{"message": {"content": json.dumps(PROPOSAL), "refusal": "No"}}]},
    {"choices": [{"message": {"content": "x" * (llm.MAX_CONTENT_BYTES + 1)}}]},
    {"choices": [{"message": {"content": json.dumps({**PROPOSAL, "explanation": KEY})}}]},
    {"choices": [{"message": {"content": json.dumps(PROPOSAL)}}], "request_echo": {"Authorization": "Bearer " + KEY}},
])
def test_envelope_failures_and_credential_echo_are_rejected(monkeypatch, envelope):
    wire(monkeypatch, raw=json.dumps(envelope).encode())
    with pytest.raises(llm.LLMSelectionError) as caught:
        select()
    assert KEY not in str(caught.value)


@pytest.mark.parametrize("status", [301, 302, 303, 307, 308, 400, 401, 403, 429, 500])
def test_provider_http_errors_never_expose_error_body_url_or_credentials(monkeypatch, status, capsys):
    failure = HTTPError(URL + "/" + KEY, status, KEY, {}, BytesIO((KEY + " private provider error").encode()))
    calls = wire(monkeypatch, error=failure)
    with pytest.raises(llm.LLMSelectionError) as caught:
        select()
    assert len(calls) == 1
    assert KEY not in str(caught.value) and "provider.example" not in str(caught.value)
    assert capsys.readouterr() == ("", "")


@pytest.mark.parametrize("error", [URLError(KEY), TimeoutError(KEY), OSError(KEY)])
def test_transport_errors_are_generic_and_never_retried(monkeypatch, error):
    calls = wire(monkeypatch, error=error)
    with pytest.raises(llm.LLMSelectionError, match="connection failed") as caught:
        select()
    assert len(calls) == 1 and KEY not in str(caught.value)


@pytest.mark.parametrize("url", ["http://provider.example/v1", "http://127.1/v1", "http://2130706433/v1",
    "http://localhost.evil/v1", "ftp://localhost/v1", "https://user:password@provider.example/v1",
    "https://provider.example/v1?api_key=" + KEY, "https://provider.example/v1#x", "https://provider.example/v1?",
    "https://provider.example/v1/../other", "https://provider.example/%2e%2e/v1", "https://provider.example\\@evil/v1",
    "https://provider.example:99999/v1", "https://provider.example:/v1", "https://provider.example/\nsecret",
    "https://[::1%25eth0]/v1"])
def test_unsafe_url_is_rejected_before_transport(monkeypatch, url):
    calls = wire(monkeypatch)
    with pytest.raises(llm.LLMSelectionError):
        select(base_url=url)
    assert calls == []


@pytest.mark.parametrize("url", ["http://localhost:11434/v1", "http://127.0.0.1:1234/v1", "http://[::1]:8000/v1",
                                "https://provider.example/custom/v1/chat/completions"])
def test_exact_local_http_or_https_endpoint_supported_without_auth(monkeypatch, url):
    endpoint = url if url.endswith("/chat/completions") else url + "/chat/completions"
    calls = wire(monkeypatch, endpoint=endpoint)
    assert select(base_url=url, api_key="") == PROPOSAL
    assert calls[0][0].get_header("Authorization") is None


@pytest.mark.parametrize("key", ["bad\r\nX-Header: value", "contains space", "bad\x00key", "x" * 4097])
def test_invalid_credentials_rejected_before_network(monkeypatch, key):
    calls = wire(monkeypatch)
    with pytest.raises(llm.LLMSelectionError):
        select(api_key=key)
    assert not calls


def test_secret_cannot_be_put_in_prompt_model_or_endpoint(monkeypatch):
    calls = wire(monkeypatch)
    for overrides in [{"model": KEY}, {"base_url": URL + "/" + KEY}]:
        with pytest.raises(llm.LLMSelectionError):
            select(**overrides)
    with pytest.raises(llm.LLMSelectionError):
        llm.select_with_llm("The key is " + KEY, LIBRARY, base_url=URL, model="test", api_key=KEY)
    assert not calls


@pytest.mark.parametrize("kwargs", [{"raw": b"x" * (llm.MAX_RESPONSE_BYTES + 1)},
                                   {"headers": {"Content-Length": str(llm.MAX_RESPONSE_BYTES + 1)}},
                                   {"raw": b"\xff"}, {"endpoint": "https://other.example/"}])
def test_response_size_encoding_and_unexpected_final_url(monkeypatch, kwargs):
    wire(monkeypatch, **kwargs)
    with pytest.raises(llm.LLMSelectionError):
        select()


def test_real_urllib_redirect_chain_is_blocked_before_second_request():
    # Exercise urllib's real handler chain without creating sockets or a server.
    calls = []
    class RedirectingTransport(HTTPHandler):
        def http_open(self, request):
            calls.append(request)
            response = Response(b"", endpoint=request.full_url,
                                headers={"Location": "http://localhost:9001/stolen"}, code=302)
            response.msg = "Found"
            response.url = request.full_url
            return response
    opener = build_opener(llm.ProxyHandler({}), llm._NoRedirect(), RedirectingTransport())
    request = Request("http://localhost:9000/v1/chat/completions", data=b"{}", headers={"Authorization": "Bearer " + KEY})
    with pytest.raises(llm.LLMSelectionError, match="redirects"):
        opener.open(request)
    assert len(calls) == 1


def test_read_deadline_is_bounded_and_no_environment_key_is_discovered(monkeypatch):
    calls = wire(monkeypatch)
    monkeypatch.setenv("OPENAI_API_KEY", KEY)
    assert select(api_key="") == PROPOSAL
    assert calls[0][0].get_header("Authorization") is None
    monkeypatch.setattr(llm.time, "monotonic", lambda: 100)
    with pytest.raises(llm.LLMSelectionError, match="timed out"):
        llm._read_response(Response(b"{}"), deadline=99)


def test_actual_description_catalog_shape_and_aliases_work(monkeypatch):
    from cof_vision.cofdraw_describe import describe_catalog
    from cof_vision.cofdraw_chemistry import TEMPLATES, library_linkages
    from cof_vision.cofdraw_styles import STYLES, PALETTES
    data = Path(__file__).resolve().parents[1] / "data/framework-drawing/catalog/current.json"
    catalog = json.loads(data.read_text(encoding="utf-8"))
    raw_library = {"monomer_presets": catalog["monomer_presets"],
        "source_recipes": [{"recipe_id": key, **value} for key, value in catalog["source_metadata"].items()],
        "linkages": library_linkages(catalog), "templates": TEMPLATES, "styles": STYLES,
        "style_controls": {"palettes": list(PALETTES)}}
    descriptor = describe_catalog(raw_library)
    proposal = {**PROPOSAL, "selection": {"source_recipe_id": "pyurea-source-layout-v1", "style_id": "source"}}
    calls = wire(monkeypatch, proposal)
    assert llm.select_with_llm("Draw Pyurea in source style", descriptor, base_url=URL, model="test-model", api_key=KEY) == proposal
    prompt_catalog = json.loads(json.loads(calls[0][0].data)["messages"][1]["content"])["catalog"]
    source = next(r for r in prompt_catalog["sources"] if r["id"] == "pyurea-source-layout-v1")
    assert "Pyurea" in source["aliases"]
    assert all(r["id"] in {key for key, row in raw_library["linkages"].items() if row["supported"]} for r in prompt_catalog["linkages"])
    assert all("smiles" not in row for row in prompt_catalog["monomers"])


def test_json_escaped_credentials_cannot_be_returned(monkeypatch):
    secret = 'quoted"key\\not-for-output'
    wire(monkeypatch, {**PROPOSAL, "explanation": secret})
    with pytest.raises(llm.LLMSelectionError) as caught:
        select(api_key=secret)
    assert secret not in str(caught.value)


@pytest.mark.parametrize("sites", [0, -1, 25, True, "2", 2.5])
def test_prompt_omits_invalid_or_unbounded_reactive_site_counts(sites):
    library = deepcopy(LIBRARY)
    library["monomers"][1]["compatible_roles"][0]["reactive_sites"] = sites
    compact, _ = llm._catalog(library)
    assert "reactive_sites" not in compact["monomers"][1]["compatible_roles"][0]
