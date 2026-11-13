# Source-table dictionary

All residue numbers use native human ETV4 coordinates. Representation metrics are dimensionless. Exact numeric precision is retained. Numerical bounds are not confidence intervals. Prompt settings, blocks and shuffle draws are not independent biological replicates.

## esm-all-fixed-blocks.tsv

30 rows.

| Column | Meaning |
|---|---|
| `model` | ESM model identifier; model spaces are not assumed commensurate |
| `block` | Post-transformer-block index (before final model-wide normalization) |
| `genotype` | Substitution genotype in native ETV4 coordinates; underscore separates substitutions |
| `M` | RMS difference between full-length and ETS+CID mutation-displacement vectors |
| `B` | RMS WT full-length versus ETS+CID displacement |
| `D` | RMS mutation-displacement vector in ETS+CID |
| `T` | M/B descriptive diagnostic |
| `M_over_D` | M/D descriptive ratio; not validated regulatory strength |
| `L` | Lower numerical sensitivity bound for M/D |
| `U` | Upper numerical sensitivity bound for M/D |
## esm-native-order-margins.tsv

6 rows.

| Column | Meaning |
|---|---|
| `model` | ESM model identifier; model spaces are not assumed commensurate |
| `block` | Post-transformer-block index (before final model-wide normalization) |
| `primary` | Whether this block was the designated primary readout |
| `native_M` | Native Y401N/Y403G M |
| `max_shuffle_M` | Largest Y401N/Y403G M over the 20 shuffled prefixes |
| `margin` | native_M minus max_shuffle_M |
| `required_8eta` | Frozen numerical comparison allowance, 8 eta |
| `passes` | Whether margin strictly exceeds required_8eta |
## esm-primary-shuffles.tsv

200 rows.

| Column | Meaning |
|---|---|
| `model` | ESM model identifier; model spaces are not assumed commensurate |
| `block` | Post-transformer-block index (before final model-wide normalization) |
| `genotype` | Substitution genotype in native ETV4 coordinates; underscore separates substitutions |
| `draw` | Synthetic-control draw identifier; not an independent biological replicate |
| `M` | RMS difference between full-length and ETS+CID mutation-displacement vectors |
## esm-primary.tsv

10 rows.

| Column | Meaning |
|---|---|
| `model` | ESM model identifier; model spaces are not assumed commensurate |
| `block` | Post-transformer-block index (before final model-wide normalization) |
| `genotype` | Substitution genotype in native ETV4 coordinates; underscore separates substitutions |
| `M` | RMS difference between full-length and ETS+CID mutation-displacement vectors |
| `B` | RMS WT full-length versus ETS+CID displacement |
| `D` | RMS mutation-displacement vector in ETS+CID |
| `T` | M/B descriptive diagnostic |
| `M_over_D` | M/D descriptive ratio; not validated regulatory strength |
| `L` | Lower numerical sensitivity bound for M/D |
| `U` | Upper numerical sensitivity bound for M/D |
## esm-wt-factorial.tsv

42 rows.

| Column | Meaning |
|---|---|
| `model` | ESM model identifier; model spaces are not assumed commensurate |
| `block` | Post-transformer-block index (before final model-wide normalization) |
| `kind` | Native/shuffle or prompt category, depending on table |
| `draw` | Synthetic-control draw identifier; not an independent biological replicate |
| `rms_J` | RMS WT factorial vector interaction contrast |
## profam-NID-shuffles.tsv

420 rows.

| Column | Meaning |
|---|---|
| `block` | Post-transformer-block index (before final model-wide normalization) |
| `prompt` | Prompt identifier linked to prompt-panel.json |
| `draw` | Synthetic-control draw identifier; not an independent biological replicate |
| `rms_q` | RMS response to adding NID under the specified prompt |
## profam-matched-norm-differences.tsv

30 rows.

| Column | Meaning |
|---|---|
| `block` | Post-transformer-block index (before final model-wide normalization) |
| `set` | Matched prompt-pair set (1–5) |
| `order` | Order of two prompting proteins: ab or ba |
| `hom_norm` | NID-response magnitude under homologous prompt |
| `nonhom_norm` | NID-response magnitude under matched nonhomologous prompt |
| `difference_of_norms` | hom_norm minus nonhom_norm; not a vector norm |
## profam-order-contrasts.tsv

45 rows.

| Column | Meaning |
|---|---|
| `block` | Post-transformer-block index (before final model-wide normalization) |
| `prompt_pair` | Pair of prompt settings compared |
| `rms_order_vector_difference` | RMS difference of NID-response vectors after prompt order reversal |
## profam-prompt-responses.tsv

93 rows.

| Column | Meaning |
|---|---|
| `block` | Post-transformer-block index (before final model-wide normalization) |
| `prompt` | Prompt identifier linked to prompt-panel.json |
| `kind` | Native/shuffle or prompt category, depending on table |
| `rms_q` | RMS response to adding NID under the specified prompt |
## profam-vector-contrasts.tsv

30 rows.

| Column | Meaning |
|---|---|
| `block` | Post-transformer-block index (before final model-wide normalization) |
| `set` | Matched prompt-pair set (1–5) |
| `order` | Order of two prompting proteins: ab or ba |
| `hom_minus_nonhom` | RMS of q_hom minus q_nonhom as vectors, not the difference of their magnitudes |
| `hom_minus_homshuffle` | RMS of q_hom minus q_homshuffle as vectors, not the difference of their magnitudes |
## published-WT-affinities.tsv

7 rows.

| Column | Meaning |
|---|---|
| `construct` | Published WT construct identifier |
| `KD_nM` | Published dissociation constant converted to nanomolar |
| `SEM_nM` | Published standard error converted to nanomolar |
| `reported_n` | Replicate count reported by the biological source |
