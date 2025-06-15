# Met-stress DEG Analysis

## 1. Overview  
Analysis of methionine-stress phenotypes in cancer cell lines to dissect homocysteine (Hcy) responses and adaptive changes in methionine-stress-resistant revertants.

---

## 2. Research Questions  
1. **Parental Hcy response**  
   Is there a common cellular response to a shift to homocysteine (Hcy) in MDA-MB-468 and 293T cells?  
2. **Maintenance in revertants**  
   Which of those shared Hcy-response pathways are preserved in the methionine-stress-resistant lines R1 (from 293T) and R8 (from 468)?  
3. **Unique Hcy induction in revertants**  
   Are there pathways uniquely induced by Hcy in the revertants compared to their parental lines—and if so, what are they?  
4. **Baseline Met differences**  
   Under methionine (Met) conditions, do R1 and R8 already exhibit altered expression profiles relative to their parents, and are those profiles similar between R1 and R8?  

---

## 3. Contrast Definitions  

| Contrast                | Condition A      | Condition B      | Purpose                                    |
|-------------------------|------------------|------------------|--------------------------------------------|
| **468_hcy_vs_met**      | 468_Hcy          | 468_Met          | Q1: Parental Hcy response in 468          |
| **293t_hcy_vs_met**     | 293T_Hcy         | 293T_Met         | Q1: Parental Hcy response in 293T         |
| **r8_hcy_vs_met**       | R8_Hcy           | R8_Met           | Q2: Maintenance of Hcy response in R8     |
| **r1_hcy_vs_met**       | R1_Hcy           | R1_Met           | Q2: Maintenance of Hcy response in R1     |
| **r8_vs_468_hcy**       | R8_Hcy           | 468_Hcy          | Q3: Unique Hcy induction in R8            |
| **r1_vs_293t_hcy**      | R1_Hcy           | 293T_Hcy         | Q3: Unique Hcy induction in R1            |
| **r8_vs_468_met**       | R8_Met           | 468_Met          | Q4: Baseline Met differences in R8        |
| **r1_vs_293t_met**      | R1_Met           | 293T_Met         | Q4: Baseline Met differences in R1        |

---

## 4. Analysis Plan

### Q1: Common Hcy response in parental lines  
- **Contrasts:**  
  - `468_hcy_vs_met`  
  - `293t_hcy_vs_met`  
- **Approach:**  
  1. Perform DE analysis for each contrast.  
  2. Identify shared up- and down-regulated genes/pathways between the two.

### Q2: Maintenance of Hcy response in revertants  
- **Contrasts:**  
  - `r8_hcy_vs_met`  
  - `r1_hcy_vs_met`  
- **Approach:**  
  1. From Q1’s shared DE gene list, test which remain DE in R8 and R1 under Hcy.  
  2. Characterize preserved pathways.

### Q3: Unique Hcy induction in revertants vs. parents  
- **Contrasts:**  
  - `r8_vs_468_hcy`  
  - `r1_vs_293t_hcy`  
- **Approach:**  
  1. Identify DE genes/pathways in each revertant vs. parent under Hcy.  
  2. Remove any genes already shared in Q1 to find unique signatures.

### Q4: Baseline Met expression differences  
- **Contrasts:**  
  - `r8_vs_468_met`  
  - `r1_vs_293t_met`  
- **Approach:**  
  1. Identify DE genes/pathways in each revertant vs. parent under Met.  
  2. Compare R8 vs. R1 DE sets to assess similarity.

---

## 5. Plotting Strategy

For each question, generate separate plots for up-regulated and down-regulated genes.

| Question(s)        | Contrasts                                                   | Upset Plot                  | Venn Diagram                 |
|--------------------|-------------------------------------------------------------|-----------------------------|------------------------------|
| **Q1 & Q2**        | `468_hcy_vs_met`, `293t_hcy_vs_met`, `r8_hcy_vs_met`, `r1_hcy_vs_met` | – Up-regulated (4 sets)<br>– Down-regulated (4 sets) | – Up-regulated<br>– Down-regulated |
| **Q3**             | `r8_vs_468_hcy`, `r1_vs_293t_hcy`                           | – Up-regulated (2 sets)<br>– Down-regulated (2 sets) | – Up-regulated<br>– Down-regulated |
| **Q4**             | `r8_vs_468_met`, `r1_vs_293t_met`                           | – Up-regulated (2 sets)<br>– Down-regulated (2 sets) | – Up-regulated<br>– Down-regulated |

- **Total plots:**  
  - 3 Upset-plots for up-regulated sets (one per question-group)  
  - 3 Upset-plots for down-regulated sets  
  - 3 Venn diagrams for up-regulated sets  
  - 3 Venn diagrams for down-regulated sets  

---

## 6. How to Run

1.  **Install Dependencies:**
    This project uses `renv` to manage R dependencies. To install the required packages, open an R session in the project root and run:
    ```R
    renv::restore()
    ```

2.  **Run the Full Pipeline:**

    You can run the pipeline in two ways:

    **A. From within an R session (recommended for interactive use):**
    
    Open an R session in the project root and run:
    ```R
    source("scripts/main.R")
    ```

    **B. From the command line (for automation):**

    Navigate to the project root in your terminal (Zsh, Bash, etc.) and execute the following command:
    ```sh
    Rscript scripts/main.R
    ```

    This will run all scripts in the correct order and generate all results in the `results/` directory and final reports in the `reports/` directory.

---

