## Relaxed-Spin Data Metadata & Streaming Format

### 1. Core Terminology & Data Layout

* 
**Node:** The logical entity represented in data lines by a pair of values: `spin_i X_i`. `spin_i` is either `+1` or `-1`, and `X_i` is a floating-point value.


* 
**Data Record (File):** A raw numeric line starting with a `time` float, followed by space-separated node pairs.


```text
0.000  1 0.12  -1 0.44  1 0.91

```


* 
**Metadata Record (File):** A machine-readable JSON object embedded in a comment line prefixed with `#@`.


* 
**Stream Message:** The exact same JSON object transmitted via NDJSON directly, omitting the `#@` file prefix.



---

### 2. File-Based Syntax & Compatibility Rationale

* 
**Legacy Parser Assumption:** Legacy parsers are assumed to ignore lines beginning with `#`. This specification deliberately embeds machine-readable metadata within these comment lines to maintain strict backward compatibility without breaking existing parsers.


* 
**Syntax Distinction:** Lines starting with `#` but not `#@` are ordinary comments. Enhanced parsers ignore ordinary comments and only interpret lines beginning with `#@`.


```text
# [cite_start]This is an ordinary comment (ignored by metadata parsers)
#@ {"type": "node", "index": 5, "name": "pass_a"}

```


* 
**Precedence:** If multiple metadata records assign the same field to the same entity, the **later record overrides the earlier one** ("last record wins").


* 
**Parser Modes:** * *Permissive Mode (Default):* Warns and continues on malformed records; ignores unknown record types and unrecognized fields.


* 
*Strict Mode:* Raises an error on any malformed metadata, unknown types, or conflicting indexing declarations.





---

### 3. Core Metadata Objects

Every metadata record requires a `"type"` field to identify its semantic meaning.

#### Format Declaration

Identifies the specific layout and versioning parameters.

```json
{"type": "format", "name": "relaxed-spins", "version": 1}

```

#### Node Indexing

Defines the global index base (`0` or `1`) for logical node references.

* Must appear before any records referencing node indices.


* It is **global** and cannot change mid-session. If omitted, parsers default to `base: 1`.


* 
**Producer Recommendation:** To avoid ambiguity between different development ecosystems (e.g., Julia vs. C++/Python toolchains), data producers **should explicitly emit** a `node_indexing` record.



**Example in a File:**

```text
#@ {"type": "node_indexing", "base": 1}

```

**Example in a Stream:**

```json
{"type": "node_indexing", "base": 1}

```

#### Node Metadata

Maps a logical node position index (not raw file columns) to semantic attributes.

```json
{"type": "node", "index": 1, "name": "CTRL_1"}

```

#### Optional Context Metadata

```json
{"type": "created", "value": "2026-01-26T18:44:00Z"}
{"type": "title", "value": "Trajectory for test run 42"}
{"type": "notes", "value": "Simulation after thermal relaxation."}

```

---

### 4. Stream-Native Formats

In streaming contexts, data samples and control mechanisms are transmitted as pure newline-delimited JSON objects without comment prefixes.

#### Stream Samples

Replaces raw text rows with structured JSON arrays.

```json
{"type": "sample", "time": 0.0, "r_spins": [{"state": [1, 0.12]}, {"state": [-1, 0.44]}, {"state": [1, 0.91]}]}

```

#### Bidirectional Control & Feedback

Enables a consumer or controller to send commands back to the data source, using an optional `"id"` for request-response tracking.

```json
{"type": "control", "id": "msg-17", "command": "set_parameter", "name": "temperature", "value": 0.25}
{"type": "ack", "id": "msg-17", "status": "ok"}
{"type": "error", "id": "msg-17", "message": "Invalid value"}

```

---

### 5. Consolidated Examples

#### Simple header with simulation and node metadata

```text
# Simulation metadata (this is an ordinary comment)
#@ {"type": "format", "name": "relaxed-spins", "version": 1}
#@ {"type": "created", "value": "2026-01-26T18:44:00Z"}
#@ {"type": "title", "value": "Trajectory for test run 42"}
#@ {"type": "notes", "value": "Simulation after thermal relaxation."}
#@ {"type": "node_indexing", "base": 1}
#@ {"type": "node", "index": 1, "name": "CTRL_1"}
#@ {"type": "node", "index": 5, "name": "pass_a"}

0.000  1 0.12  -1 0.44 ...
0.100  1 0.15  -1 0.43 ...

```

#### Complete Minimal Stream (NDJSON)

```json
{"type": "format", "name": "relaxed-spins", "version": 1}
{"type": "node_indexing", "base": 1}
{"type": "node", "index": 1, "name": "CTRL_1"}
{"type": "sample", "time": 0.0, "node_states": [[1, 0.12], [-1, 0.44]]}

```