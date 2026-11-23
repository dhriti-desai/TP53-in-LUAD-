# imports panads library and uses short form <<pd>> so easier to call
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter

# data set 1: data_mutations

# \t tells pd to seperate by tabs
mutations = pd.read_csv("data_mutations.txt", sep="\t")
mutations.to_csv("mutations.csv", index=False)

"""
gene_counts = mutations["Gene"].value_counts()

top_10 = gene_counts.head(10)
top_10.plot(kind="bar")

plt.title("Top 10 Most Frequenlty Mutates Genes")
plt.xlabel("Gene")
plt.ylabel("Number of Mutations")

plt.show()
"""

"""
# graph showing distribution on mutation types
dist_mut = mutations["Variant_Classification"].value_counts()
dist_mut_10 = dist_mut.head(10)
dist_mut_10.plot(kind="bar")

plt.title("Mutation Distribtuion")
plt.xlabel("Mutation Type")
plt.ylabel("Frequency")
plt.show()
"""

"""
print(mutations["Tumor_Sample_Barcode"].head())
print(mutations["Tumor_Sample_Barcode"].nunique())
"""

grouped = mutations.groupby("Hugo_Symbol")
gene_frequency = grouped["Tumor_Sample_Barcode"].nunique()
total_samples = mutations["Tumor_Sample_Barcode"].nunique()

"""
# sort_values() sorts them from high to low --> and asdending=False means it will be from high to low not low to high
top_genes = gene_frequency.sort_values(ascending=False)
top_10_gene = top_genes.head(10)
top_10_gene.plot(kind="bar")

plt.title("Top 10 Genes")
plt.xlabel("Gene")
plt.ylabel("Number of Tumour Samples")
plt.show()
"""

prevalence = (gene_frequency/total_samples)*100

top_10_prevalence = prevalence.sort_values(ascending=False).head(10)
top_10_prevalence.to_csv("top_10_gene_prevalence.csv")

# data set 2 --> clinical_patient
clinical = pd.read_csv("data_clinical_patient.txt", sep="\t", skiprows=4)

"""
print(clinical["OS_STATUS"].value_counts(dropna=False)) # number of patients in survival category
print(clinical["OS_MONTHS"].describe()) # gives count, mean, standard deviation, minimum, and quartiles of survival time
print(clinical["OS_MONTHS"].isna().sum()) # counts missing survival-time values
"""

# dropna() removes rows containing missing values --> does not direclty change the csv only if told to
survival = clinical.dropna(subset=["OS_STATUS", "OS_MONTHS"])

"""
print(survival.shape)
print(survival["OS_STATUS"].value_counts())
"""
# 182 deaths observed and 323 living --> total = 505
"""
survival["event"] creates new column in data set
since the data is 1:LIVING or 0: DECEASED; map removes the words and js keeps 0 or 1

death categories: 0 = censored | 1 = observed
"""
survival["event"] = survival["OS_STATUS"].map({"0:LIVING": 0, "1:DECEASED": 1})

"""
print(survival[["OS_STATUS", "OS_MONTHS", "event"]].head())
print(survival["event"].value_counts())
"""
kmf = KaplanMeierFitter()

"""
kmf.fit(durations=survival["OS_MONTHS"], event_observed=survival["event"], label="All LUAD Patients")

kmf.plot_survival_function()

plt.title("Kaplan-Meier Survival Curve")
plt.xlabel("Time (Months)")
plt.ylabel("Estimated Survival Probability")
plt.show()
"""

clin_samp = pd.read_csv("data_clinical_sample.txt", sep="\t", skiprows=4)

"""
checking if each gene corresponds to patient id --> they all corresspond

# get all unique sample IDs from the mutation dataset
mutation_samples = mutations["Tumor_Sample_Barcode"].dropna().unique()

# get all unique sample IDs from the clinical sample dataset
clinical_samples = clin_samp["SAMPLE_ID"].dropna().unique()

# count how many mutation sample IDs appear in the clinical dataset
matching_samples = pd.Series(mutation_samples).isin(clinical_samples)

print("Mutation samples:", len(mutation_samples))
print("Matching samples:", matching_samples.sum())
print("Unmatched samples:", (~matching_samples).sum())
"""

mutation_patients = mutations.merge(clin_samp[["SAMPLE_ID", "PATIENT_ID"]], left_on="Tumor_Sample_Barcode", right_on="SAMPLE_ID", how="left")
mutation_survival = mutation_patients.merge(survival[["PATIENT_ID", "OS_MONTHS", "OS_STATUS", "event"]], on="PATIENT_ID", how="inner")

# keep only patient ID and gene
patient_gene = mutation_survival[["PATIENT_ID", "Hugo_Symbol"]]

# remove duplicate patient-gene combinations
patient_gene = patient_gene.drop_duplicates()

# count distinct patients for each gene
patient_gene_counts = (patient_gene.groupby("Hugo_Symbol")["PATIENT_ID"].nunique().sort_values(ascending=False))
top_10_patients = patient_gene_counts.head(10)

# filter mutation_survival to rows where Hugo_Symbol equals "TP53"
tp53_mutations = mutation_survival[mutation_survival["Hugo_Symbol"] == "TP53"]

tp53_patient_count = tp53_mutations["PATIENT_ID"].nunique() # count unique patient IDs

tp53_patients = mutation_survival[ mutation_survival["Hugo_Symbol"] == "TP53"]["PATIENT_ID"].unique()
survival["TP53_mutated"] = survival["PATIENT_ID"].isin(tp53_patients).astype(int)


# Does the estimated survival probability differ between patients with and without a recorded TP53 mutation?
mutated = survival[survival["TP53_mutated"] == 1]
not_mutated = survival[survival["TP53_mutated"] == 0]

"""
kmf.fit(mutated["OS_MONTHS"],event_observed=mutated["event"],label="TP53 Mutated")
ax = kmf.plot_survival_function()

kmf.fit(not_mutated["OS_MONTHS"], event_observed=not_mutated["event"], label="No Recorded TP53 Mutation")
kmf.plot_survival_function(ax=ax)

plt.title("Overall Survival by TP53 Mutation Status")
plt.xlabel("Time (Months)")
plt.ylabel("Estimated Survival Probability")
plt.show()
"""

# missing patients
survival_patients = set(survival["PATIENT_ID"])
mutation_patient_ids = set(mutation_survival["PATIENT_ID"])
missing_patients = survival_patients - mutation_patient_ids

survival_profiled = survival[ ~survival["PATIENT_ID"].isin(missing_patients)].copy()
mutated = survival_profiled[survival_profiled["TP53_mutated"] == 1]
not_mutated = survival_profiled[survival_profiled["TP53_mutated"] == 0]

# log rank test
results = logrank_test(mutated["OS_MONTHS"], not_mutated["OS_MONTHS"], event_observed_A=mutated["event"], event_observed_B=not_mutated["event"])

# cox proportional hazards regression
cox_data = survival_profiled[["OS_MONTHS", "event", "TP53_mutated"]].copy()
cph = CoxPHFitter()
cph.fit(cox_data, duration_col="OS_MONTHS", event_col="event")

# checks if cph assumptions are reasonable 
cph.check_assumptions(cox_data, p_value_threshold=0.05, show_plots=True)


#checking if missing patients changes our data
# log test
cox_data_all = survival[["OS_MONTHS", "event", "TP53_mutated"]].copy()
cph_all = CoxPHFitter()

cph_all.fit(cox_data_all, duration_col="OS_MONTHS",event_col="event")

print(cph_all.summary)

# cox test
mutated_all = survival[survival["TP53_mutated"] == 1]
not_mutated_all = survival[survival["TP53_mutated"] == 0]

results_all = logrank_test(mutated_all["OS_MONTHS"], not_mutated_all["OS_MONTHS"], event_observed_A=mutated_all["event"], event_observed_B=not_mutated_all["event"])

survival_adjusted = survival.copy()

# keep patients with complete age and stage information
survival_adjusted = survival_adjusted.dropna(subset=["AGE", "AJCC_PATHOLOGIC_TUMOR_STAGE"]).copy()

def group_stage(stage):
    if stage.startswith("STAGE III"):
        return "III"
    elif stage.startswith("STAGE IV"):
        return "IV"
    elif stage.startswith("STAGE II"):
        return "II"
    elif stage.startswith("STAGE I"):
        return "I"

survival_adjusted["STAGE_GROUP"] = (
    survival_adjusted["AJCC_PATHOLOGIC_TUMOR_STAGE"].apply(group_stage)
)

stage_dummies = pd.get_dummies(survival_adjusted["STAGE_GROUP"], prefix="STAGE", drop_first=True, dtype=int)

cox_data_adjusted = pd.concat([survival_adjusted[["OS_MONTHS", "event", "TP53_mutated", "AGE"]], stage_dummies], axis=1)

cph_adjusted = CoxPHFitter()

cph_adjusted.fit(cox_data_adjusted, duration_col="OS_MONTHS", event_col="event")

# cph_adjusted.check_assumptions(cox_data_adjusted, p_value_threshold=0.05, show_plots=True)

# two hypothetical patient profiles --> compares with age + tumour sample vs. without against survival rates
profiles = pd.DataFrame({"TP53_mutated": [0, 1], "AGE": [survival_adjusted["AGE"].median()] * 2, "STAGE_II": [0, 0], "STAGE_III": [0, 0], "STAGE_IV": [0, 0]})

# predict survival curves
adjusted_curves = cph_adjusted.predict_survival_function(profiles)

# only makes one graph
fig, ax = plt.subplots(figsize=(10, 6))

# plot both curves on the same axes
ax.plot(adjusted_curves.index, adjusted_curves.iloc[:, 0], label="No recorded TP53 mutation")

ax.plot(adjusted_curves.index, adjusted_curves.iloc[:, 1],label="TP53 mutated")

ax.set_title("Age- and Stage-Adjusted Survival by TP53 Status", fontsize=14)
ax.set_xlabel("Time (Months)")
ax.set_ylabel("Predicted Survival Probability")
ax.set_ylim(0, 1)
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()