# Project documents not in this site

`docs/` also holds working documents that are not part of the user guide — design scopes,
campaign plans, and the engineering audit. They are listed here so nothing in the
repository is unreachable from the site, but they are not rendered as chapters: several
are point-in-time records rather than current instructions, and folding them into the
nav would present them as guidance.

| Document | What it is |
|---|---|
| [`claude_engineering_audit.md`](https://github.com/ncsa/eidolon/blob/develop/docs/claude_engineering_audit.md) | Defect taxonomy and case histories behind the project's vetting rules |
| [`access_report_draft.md`](https://github.com/ncsa/eidolon/blob/develop/docs/access_report_draft.md) | ACCESS allocation report draft — benchmark and validation numbers |
| [`sv_polish_roadmap.md`](https://github.com/ncsa/eidolon/blob/develop/docs/sv_polish_roadmap.md) | Structural-variant work plan |
| [`scn_status.md`](https://github.com/ncsa/eidolon/blob/develop/docs/scn_status.md), [`scn_phase2_af_design.md`](https://github.com/ncsa/eidolon/blob/develop/docs/scn_phase2_af_design.md) | Subclonal copy-number status and allele-fraction design |
| [`transition_matrix_validation_plan.md`](https://github.com/ncsa/eidolon/blob/develop/docs/transition_matrix_validation_plan.md) | Validation plan for the SNP transition matrix |
| [`subcontig_chunking_plan.md`](https://github.com/ncsa/eidolon/blob/develop/docs/subcontig_chunking_plan.md) | Design notes for sub-contig chunking |
| [`longread_epic_scope.md`](https://github.com/ncsa/eidolon/blob/develop/docs/longread_epic_scope.md) | Scope for long-read support ([#319](https://github.com/ncsa/eidolon/issues/319)) |
| [`rename_eidolon_scope.md`](https://github.com/ncsa/eidolon/blob/develop/docs/rename_eidolon_scope.md) | Scope of the `rneat` → `eidolon` rename |

Model provenance — what each shipped default is, where it came from, and whether it has
been measured — lives in
[`eidolon-core/src/models/model_data/README.md`](https://github.com/ncsa/eidolon/blob/develop/eidolon-core/src/models/model_data/README.md).
