# JSON Data Specification: Stabilized Record Types

## 1. Overview

This specification documents the stabilized `"type"` fields used for storing and exchanging Dice-related data. Every metadata and data record in this specification is represented as a JSON object containing a discriminator field `"type"`, which identifies the object's schema and semantic role.

---

## 2. Base Requirements

Every valid JSON record in this specification MUST conform to the following base requirements:

* **Type Discriminator:** A top-level string property named `"type"` is required in every JSON object.
* **Parsing Policy:** In Permissive Mode, parsers ignore unrecognized fields within a known record type. In Strict Mode, unrecognized or malformed fields cause parsing failures.

---

## 3. Stabilized Record Types

### 3.1. Format Declaration (`"type": "format"`)

Identifies the specific layout format and layout version of the dataset.

* **Type Value:** `"format"`
* **Required Attributes:**
  * `name` (string): The identifier of the data format layout specification (e.g., `"relaxed-spins"`).
  * `version` (integer): The version number of the format layout (e.g., `1`).

**JSON Example:**
```json
{"type": "format", "name": "relaxed-spins", "version": 1}
```

---

### 3.2. Node Indexing (`"type": "node_indexing"`)

Defines the global index base (`0` or `1`) for logical node references across the session or dataset.

* **Type Value:** `"node_indexing"`
* **Required Attributes:**
  * `base` (integer): Global index base for node positions. Must be `0` (0-based indexing) or `1` (1-based indexing).
* **Behavioral Rules:**
  * Must appear before any metadata or sample records that reference node indices.
  * The index base is **global** and immutable mid-session.
  * If omitted, parsers default to `base: 1`. Data producers are strongly recommended to emit this record explicitly to avoid multi-language ambiguity.

**JSON Example:**
```json
{"type": "node_indexing", "base": 1}
```

---

### 3.3. Node Metadata (`"type": "node"`)

Maps a logical node position index to semantic attributes such as human-readable names.

* **Type Value:** `"node"`
* **Required Attributes:**
  * `index` (integer): The logical node position index, evaluated against the active `node_indexing` base.
  * `name` (string): Semantic identifier or string label assigned to the node (e.g., `"CTRL_1"`).

**JSON Example:**
```json
{"type": "node", "index": 1, "name": "CTRL_1"}
```

---

### 3.4. Creation Timestamp (`"type": "created"`)

Provides contextual metadata regarding when the simulation run or dataset was generated.

* **Type Value:** `"created"`
* **Required Attributes:**
  * `value` (string): Timestamp string (typically ISO 8601 formatted, e.g., `"2026-01-26T18:44:00Z"`).

**JSON Example:**
```json
{"type": "created", "value": "2026-01-26T18:44:00Z"}
```

---

### 3.5. Dataset Title (`"type": "title"`)

Provides a human-readable title or header for the dataset or simulation run.

* **Type Value:** `"title"`
* **Required Attributes:**
  * `value` (string): Title string describing the dataset (e.g., `"Trajectory for test run 42"`).

**JSON Example:**
```json
{"type": "title", "value": "Trajectory for test run 42"}
```

---

### 3.6. Annotations and Notes (`"type": "notes"`)

Contains detailed annotations, operational notes, or execution comments.

* **Type Value:** `"notes"`
* **Required Attributes:**
  * `value` (string): Free-form text note (e.g., `"Simulation after thermal relaxation."`).

**JSON Example:**
```json
{"type": "notes", "value": "Simulation after thermal relaxation."}
```

---

### 3.7. Stream Sample (`"type": "sample"`)

Represents a relaxed spin configuration.

* **Type Value:** `"sample"`
* **Required Attributes:**
  * `time` (number): The timestamp or simulation time offset (floating-point or integer).
  * `r_spins` (array of objects): Array of spin state entries. Each element contains a `state` object property formatted as a two-element array `[spin_i, X_i]`, where `spin_i` is `1` or `-1` and `X_i` is a floating-point value.

**JSON Example:**
```json
{"type": "sample", "time": 0.0, "r_spins": [{"state": [1, 0.12]}, {"state": [-1, 0.44]}, {"state": [1, 0.91]}]}
```

---

## 4. Stabilized Schema Summary

| `type` Value | Category | Fields & Types | Description |
| :--- | :--- | :--- | :--- |
| `format` | Specification Header | `name` (string), `version` (int) | Identifies data format name and version number. |
| `node_indexing` | Data Alignment | `base` (int: 0 or 1) | Sets global 0-based or 1-based indexing rule. |
| `node` | Entity Attributes | `index` (int), `name` (string) | Binds a logical node index to a named attribute. |
| `created` | Context | `value` (string) | ISO 8601 creation timestamp. |
| `title` | Context | `value` (string) | Dataset or trajectory title. |
| `notes` | Context | `value` (string) | Explanatory notes or annotations. |
| `sample` | Data Payload | `time` (number), `r_spins` (array) | Structured data sample object for streaming. |

---

## 5. Usage in Serialization Contexts

```json
{"type": "format", "name": "relaxed-spins", "version": 1}
{"type": "node_indexing", "base": 1}
{"type": "node", "index": 1, "name": "CTRL_1"}
{"type": "sample", "time": 0.0, "r_spins": [{"state": [1, 0.12]}, {"state": [-1, 0.44]}, {"state": [1, 0.91]}]}
```
