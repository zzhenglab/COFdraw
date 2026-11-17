# Describe a COF with rules or an optional model API

COFdraw can turn a description into a catalog-bound drawing plan. Rules mode runs locally and is the default. Optional LLM mode asks a user-configured provider to select existing IDs, then sends those IDs through the same deterministic chemistry planner. Both modes return a plan for review; interpreting a description does not render, create a crystal, or execute model-generated chemistry.

## Local application contract

`GET /api/drawing/describe/catalog` returns examples and selectable source, monomer, linkage, template, style and palette IDs. It contains descriptive labels/aliases and qualifications, without executable recipes or monomer SMILES. `GET /api/drawing/science/capabilities` exposes the current description request schema.

`POST /api/drawing/describe` accepts:

```json
{
  "text": "Draw Pyurea with building-block highlights and a soft palette",
  "mode": "rules",
  "selection": null,
  "llm": null
}
```

Descriptions are limited to 4,000 characters. The response is the deterministic planner's status, selected parameters, evidence, assumptions and unresolved items. A ready source plan identifies its registered recipe; a ready design supplies a bounded `/prepare` request. Review that plan before using the existing source/prepare/render workflow. Recognized unsupported chemistry, conflicting constraints and unresolved items prevent a ready plan.

Manual choices can be submitted in `selection`, using only `source_recipe_id`, `monomer_a_id`, `monomer_b_id`, `linkage`, `template_id`, `style_id` and `palette`. A source ID cannot be combined with new design chemistry. Revalidate manual choices with `mode: "rules"`; this requires no provider call. With explicit manual or LLM-selected IDs, unmatched residual words become review warnings rather than automatically blocking an otherwise valid plan. Recognized chemistry/geometry conflicts remain blocking. Review these warnings: valid catalog parameters do not guarantee that every phrase was interpreted correctly.

For optional API interpretation, use `mode: "llm"`, leave `selection` null, and supply `llm: {"base_url": "https://provider.example/v1", "model": "your-chosen-model", "api_key": "..."}` in the explicit request. Enter the actual credential only in the desktop API field or an in-memory local request. The application returns a safe interpretation summary and proposal; unresolved provider items prevent drawing readiness. The provider does not override the planner's constraints.

## Provider protocol

The adapter sends one nonstreaming HTTP POST to `<base_url>/chat/completions`. A base URL already ending in `/chat/completions` is also accepted. OpenAI's base URL is `https://api.openai.com/v1`; local compatible servers often expose their own `/v1` prefix. Supply the exact prefix and a model accepted by that provider. The request contains `model`, `messages` and `stream: false`. It uses ordinary system/user text messages and does not assume support for strict JSON response modes, tool calling, temperature or model-specific generation parameters. The response must provide one text completion in `choices[0].message.content`. [Chat Completions reference](https://developers.openai.com/api/reference/resources/chat/subresources/completions/methods/create).

The instructions request one JSON object:

```json
{
  "selection": {"source_recipe_id": "pyurea-source-layout-v1", "style_id": "source"},
  "explanation": "The named source matches the requested drawing.",
  "unresolved": []
}
```

Only exact catalog IDs are accepted. Unsupported linkages, unknown IDs, source/design mixtures, extra keys, null IDs, fenced output, malformed/duplicate-key JSON, tool calls and truncated completions are rejected. Explanations are limited to 600 characters and unresolved items to eight strings of 300 characters each. Partial valid choices can still require clarification from the deterministic planner.

Model/provider parameter support varies; a provider advertising compatibility must still implement this text Chat Completions shape. The adapter reports bounded generic errors instead of displaying provider response bodies. Automated tests use mocked provider transport, including urllib's real redirect handler chain. No live API key or paid-provider test is included in the repository.

## Credential and data handling

Nonempty keys travel only in the HTTP `Authorization: Bearer ...` header, as specified by the [API authentication reference](https://developers.openai.com/api/reference/overview#authentication). An empty key supports local providers without authentication. The adapter does not discover environment credentials, persist configuration, log requests or retry automatically. Keys containing whitespace/control characters are rejected; a key appearing in request content or a provider response also causes rejection.

Only the description and whitelisted catalog IDs, labels, aliases, role names/counts, finite-template choices and source qualifiers enter the prompt. Saved projects, prior reports, raw recipes, monomer SMILES, filesystem paths and credentials are not collected. The description itself is sent to the selected provider, whose policies and account usage apply. For first-time OpenAI account/API setup, use the [official quickstart](https://developers.openai.com/api/docs/quickstart); COFdraw takes explicitly supplied credentials rather than reading ambient environment settings.

Remote endpoints require HTTPS. HTTP is accepted only with the exact hosts `localhost`, `127.0.0.1` or `::1`. URL credentials, query strings, fragments, encoded/ambiguous path forms and control characters are rejected. Redirects are blocked before a second request; authentication is never forwarded to a redirected destination. Environment proxies are disabled. TLS verification uses Python's normal HTTPS handling.

Requests are limited to 192 KiB, provider responses to 64 KiB and completion text to 8 KiB. The socket timeout is 30 seconds and response reads use a 45-second deadline with the remaining socket timeout reduced as needed. Failed requests are not retried. Provider errors return safe messages without the endpoint, key or response body. Description plans remain transient in the local service and do not enter its scientific report history.

## Agent helper

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/description-options describe-catalog
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/description-plan describe --request description.json
```

The helper saves the returned plan in `result.json` and performs no drawing. It does not copy the request file, API settings or key into the output directory, manifests or logs. A request file created by the caller remains the caller's file; prefer the desktop field or an in-memory local request for credentials, and keep credential-bearing files out of source control and shared output folders. After reviewing a ready plan, use the existing drawing helper commands to produce actual exports.
