# AI career anxiety scale: R code

Analysis code for the manuscript *Measuring AI career anxiety in university students: initial validation of a nine-item scale in Vietnam* (under review).

## Data

Individual-level data are not posted. Participants were told that results would be reported only in aggregate, so de-identified data (the 464 cleaned responses used in every analysis) may be shared on reasonable request only if the research ethics committee that reviewed the study confirms that this is consistent with the consent given. `data/codebook.csv` describes each variable. Item wording in Vietnamese and English is given in Table S9 of the Supplementary Material.

Responses were collected from 27 to 30 October 2026, after ethics approval, through an online Google Forms survey of undergraduates at three universities in Vietnam. Of 583 responses received, 6 duplicate submissions (latest kept), 82 failed attention checks and 31 straightlined responses were removed (`ket_qua/bang/Table0_participant_flow.csv`). The form did not record each student's university.

The raw export held e-mail addresses that were used only to block duplicate submissions; it is no longer kept, and `00_Lam_sach.R` documents the cleaning steps that were applied to it. Because of this, the attention-check sensitivity analysis (Table S11, script 08) uses a response sheet without e-mail addresses, in which duplicates cannot be identified: it contains 83 failed responses against 82 in the participant flow.

## Ethics

Ethics approval was granted by the Research Ethics Committee for Social and Behavioural Sciences, Institute for Educational Research and Knowledge Transfer, Vietnam (Approval No. 02/HĐ-REK.IRB; Protocol No. REKIRB-0226; 27 October 2026).

## Reproducing the results

With the data file placed at `data/ai_career_anxiety_oct2026_n464.csv`, open `CHAY_TAT_CA.R` in RStudio and click Source, or run `Rscript CHAY_TAT_CA.R` from this folder. Required packages: lavaan, semTools, psych, jsonlite, ggplot2, patchwork and ragg (readxl only for the raw-data cleaning step). Tables are written to `ket_qua/bang/`, figures to `ket_qua/hinh/` and `ket_qua/hinh_journal/`, and full text output to `ket_qua/*.txt`.

| Script | Content |
|---|---|
| `00_CONFIG.R` | File paths and settings |
| `00_Lam_sach.R` | Cleaning of the raw export (runs only if the raw file is present) |
| `01_Phan_tich.R` | Descriptives, EFA (half A), CFA (full sample and half B), reliability, discriminant validity, invariance |
| `02_Figures.R`, `03_Figures_journal.R` | Figures |
| `04_Phan_tich_bo_sung.R` | Omega, WLSMV models, latent means, MIMIC |
| `05_Phan_tich_chat_luong.R` | Differential relations, IRT conversion, norms |
| `06_Phan_tich_phan_bien.R` | Polychoric EFA, ESEM, local dependence, ordinal invariance, MIMIC confidence intervals |
| `08_S11_attention_check_from_sheet.R` | Sensitivity analysis for respondents who failed the instructed-response item (Table S11; needs the response sheet) |
| `07_Kiem_tra_truoc_khi_cong_bo.R` | Pre-release audit: data structure, absence of personal data, participant flow, and independent recomputation of key reported values |

Random splits use `set.seed(2026)`. Code comments are in Vietnamese.

## Licence

Code: MIT (see `LICENSE`).
