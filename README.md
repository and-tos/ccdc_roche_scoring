# ccdc_roche_scoring

Contains code for template based docking and RF-PLP scoring function.

Please cite:

Tosstorff, A. et al. A high quality, industrial data set for binding affinity prediction: performance comparison in different early drug discovery scenarios. J. Comput. Aid. Mol. Des. 36, 753–765 (2022).

---
**Installation**
```

conda env create -f environment.yaml -p /path/to/env/ccdc_roche_scoring
conda activate ccdc_roche_scoring
poetry install

```

**Pre-commit hook**

In order to get in line with current commit conventions of the team, developers can set up a **pre-commit** hook that is triggered before every commit.

To first install it:

```bash
# on the conda env
$> pip install pre-commit
$> pre-commit install
...
$> git commit -am "message" # the hooks will trigger
```

---
**Template based docking in Python**

```
# Folder Structure should look like this:
docking_dir
 |-apo_proteins: contains protein and water files: strucid_dry.mol2, strucid_dry.pdb, strucid_water_0.mol2,...
 |-tmp_aligned_3d_sdf_sanitized
 | # Contains the ligand templates: native_ligands.sdf
 | # Manual curation of the templates is highly recommended to ensure proper tautomers, avoid low quality structures
 |
 |-docking: contains docking_job_0,...

# Dock ligands.sdf, with default target settings. Treat Met713 and Met712 flexibly
cd docking_dir/docking
docking.py --input_ligands ligands.sdf -t default -fr=Met713 Met712

# Concatenate results
join_docked_rf_counts.py -t [target]
```

---
**Template based docking from MOE**

One can either dock into MOE, one of the prepared projects or into an MDB of binding sites.
The binding sites MDB should have the columns
- 'moe': receptor-ligand complex
- 'ligand': 2D ligand
- 'code': STRUCID

If one of the defined projects is selected, the return results will include the RF-PLP score.

---

**Generating project specific parameters for RF-PLP:**

- First, run template docking from MOE checking the verbose box. Docked ligands have to contain attribute with pIC50.
- In the docking directory run stat_potential.py

The script will check for the presence of series definitions. If '../../series_assignment.csv' exists, series
will be assigned and parameters will be fitted per series. Alternatively, a CSV with columns SRN and split can
be passed. The model will then be fitted to the train and validation compounds, not the test compounds.

```
# Change to your docking directory
cd /path/to/docking_dir
mkdir first iteration

# Collect results from template docking
join_docked_rf_counts.py

# fit RF-PLP
cd first_iteration
stat_potential.py --target [target_name] --task [task_name_pic50]

```

---

**Score with RF-PLP:**

- Series definitions should be stored in ccdc_roche_scoring/series_definitions/[target].json
- Model parameters should be stored in ccdc_roche_scoring/scoring_parameters/[target]

```
ccdc_roche_scoring/stat_potential_inference.py -l ligand.sdf -p protein.pdb -t target_name
```

---

**Applications for MOE**
- Start a webserver on the sHPCS:
```
moeweb -load soap_scoring_client.svl -load soap_template_docking_client.svl
```
- Adapt the SVL scripts to point to the webserver: \
`const SERVER_URL = 'server address';`


---

**Scoring in MOE**

Score an MDB of ligand poses for the active protein.

---
