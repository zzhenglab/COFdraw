# Drawing conversations

The conversation is a session-only way to review catalog drawing choices. Local rules are the default. Each new drawing or appearance change needs a ready plan and an explicit **Create drawing** or **Apply appearance** action. A successful render is reported only after the actual renderer returns exports.

Try a sequence such as:

1. `Draw TFB and TAPB`
2. `Make it pastel`
3. `Change B to DETH`

The final step rechecks real reactive groups and selects a compatible implemented linkage/template. You can also provide one catalog monomer at a time, replace an explicitly named current block, change a supported style, or adjust the exact choices in the review card. A new chemistry request creates a new design; it does not inherit a source material identity or manual coordinates.

**Use current drawing** is different: an appearance-only request returns a local intent bound to that drawing's revision. The browser preserves the existing recipe and coordinates, checks the binding again when the response arrives and when Apply is clicked, and renders that recipe directly. This works for custom/imported finite drawings without sending their SMILES, graph or coordinates to the chat endpoint or an optional provider. A changed current drawing invalidates the pending appearance plan.

Unknown names, unsupported chemistry, negation and physical property targets remain unresolved. Recognized functional counts and instance counts carry into subsequent checks. A later color request cannot satisfy an earlier physical pore-size or crystal request. **Start new** explicitly clears conversation constraints and history while keeping the canvas drawing. Local rules support bounded catalog phrases, not unrestricted scientific interpretation.

Questions such as `How do I check SA scores?`, `How do I measure pores?`, or `Where can I make 3D layers?` give qualified local guidance and buttons that open the existing tools. They do not calculate or invent a value. Scientific tools still require their own explicit structures, inputs and checks.

## API and privacy

`POST /api/drawing/describe/chat` accepts:

```json
{
  "text": "Make it pastel",
  "mode": "rules",
  "context": {
    "origin": "plan",
    "selection": {"monomer_a_id": "tfb", "monomer_b_id": "tapb"},
    "pending": [],
    "requirements": []
  }
}
```

It returns the existing plan fields plus `operation` (`drawing`, `style_only`, `help`, or `reset`) and `conversation.context` for the next turn. `style_only` has no `prepare_request` or source-reload request; its `context_token` must match the caller's local drawing revision. `help` returns allowlisted action IDs, not executable links. The original `/describe` endpoint is unchanged.

Context accepts only catalog IDs, bounded unresolved messages, recognized count/role phrases, an origin and an opaque revision token. Its schema rejects recipe, SMILES, coordinate, file and arbitrary project fields. Recognizable pasted CIF/project JSON is rejected before a provider call; use the dedicated import controls instead. The browser retains at most eight prior turn summaries in memory; it does not save chat history in project files, local storage or session storage. Reloading the application clears it.

Optional model mode must be selected and sent explicitly. It uses the existing bounded OpenAI-compatible selector and passes the text the user submits, catalog descriptors and allowlisted selected IDs. It does not attach previous conversation text, unresolved local notes, the revision token, saved projects, molecular coordinates, CIFs or local files. Credentials stay in the password field and request authorization header; validation errors are redacted, and credentials are excluded from history and exports. Selecting settings or editing a message does not contact a provider. Local help does not spend a provider request even when model mode is selected.

A provider can propose existing IDs only. Local validation checks chemical preparation and recognized explicit constraints. Wording interpreted beyond local rules remains a visible review warning; this is not a guarantee of semantic fidelity or scientific prediction.

## Building-block order

The manual picker's A and B represent the user's selected order. Selecting B no longer moves or replaces A when the constructor needs the opposite reactive-role order. The frontend checks compatible orientations and maps SMILES/labels into the required preparation roles at submission time. Inspectors show the actual role for each visible block. Library previews remain ordinary 2D chemical depictions without CSS perspective transforms.
