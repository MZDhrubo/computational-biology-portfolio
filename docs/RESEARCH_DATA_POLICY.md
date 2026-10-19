# Public research-data policy

Before every public commit, review the staged diff and confirm that each file is suitable for release.

Do not commit:

- unpublished sequence datasets, candidate lists, alignments, phylogenetic trees, or research outputs;
- raw sequencing data or private sample identifiers;
- unpublished docking tables, compound identities, model files, or simulation trajectories;
- collaborator or employer files that are not cleared for public release;
- credentials, API tokens, passwords, personal local paths, or environment files.

Use synthetic examples for public demonstrations. `.gitignore` reduces accidental staging but is not a security boundary. If sensitive information is accidentally committed, removing it in a later commit may not remove it from repository history; rotate exposed credentials and follow GitHub's sensitive-data removal guidance where needed.
