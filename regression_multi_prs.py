from datetime import datetime
from pathlib import Path
import subprocess
from configparser import ConfigParser


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve
import statsmodels.api as sm
import os
from scipy.stats import mannwhitneyu, ttest_ind
from scipy.stats import gaussian_kde


script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, "config.ini")
config = ConfigParser()
config.read(config_path)

project_name = config.get("Project", "PROJECT_NAME")
output_dir = config.get("Paths", "OUTPUT_DIR")
sumstats = config.get("Paths", "SUMSTATS_FILE_ANCESTRY").split(",")


def is_none_value(value):
    if value is None:
        return True
    return str(value).strip().lower() == "none"


pheno_file = config.get("PRS", "PHENO_FILE")
covariate_file = config.get("PRS", "COVARIATE_FILE")
covariate_to_include_raw = config.get("PRS", "COVARIATE_TO_INCLUDE")

covariates_to_include = []
if not is_none_value(covariate_to_include_raw):
    covariates_to_include = [c.strip() for c in str(covariate_to_include_raw).split(",") if c.strip()]

combined_output_dir = os.path.join(output_dir, "combined")
os.makedirs(combined_output_dir, exist_ok=True)


def write_summary_file(output_path, lines):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def merge_by_common_ids(left, right):
    left = left.copy()
    right = right.copy()

    for col in ["FID", "IID"]:
        if col in left.columns:
            left[col] = left[col].astype(str)
        if col in right.columns:
            right[col] = right[col].astype(str)

    if {"FID", "IID"}.issubset(left.columns) and {"FID", "IID"}.issubset(right.columns):
        merged = left.merge(right, on=["FID", "IID"], how="inner")
        if not merged.empty:
            return merged

    if "IID" in left.columns and "IID" in right.columns:
        return left.merge(right, on="IID", how="inner")

    raise ValueError("Não foi possível identificar colunas comuns de identificação para o merge.")


dfs = []

for info_sumstat in sumstats:
    ancestry_name = info_sumstat.split(":")[-1]
    prs_best = f"{output_dir}/{ancestry_name}/{project_name}_{ancestry_name}_PRS_run.best.normalized"
    df = pd.read_csv(prs_best, sep=r"\s+")

    df = df.rename(columns={"Normalized_PRS": f"PRS_{ancestry_name}"})
    df = df[["FID", "IID", f"PRS_{ancestry_name}"]]
    dfs.append(df)


prs_df = dfs[0]

for df in dfs[1:]:
    prs_df = prs_df.merge(df, on=["FID", "IID"], how="outer")

prs_df.to_csv(f"{combined_output_dir}/{project_name}_combined_PRS_all_ancestries.txt", sep="\t", index=False)

prs_cols = [col for col in prs_df.columns if col.startswith("PRS_")]

#ate aqui eh so para preparar o arquivo
print(prs_cols)


pheno = pd.read_csv(pheno_file, sep=None, engine="python")

analysis_df = merge_by_common_ids(prs_df, pheno)

if not is_none_value(covariate_file):
    covariates = pd.read_csv(covariate_file, sep=None, engine="python")
    analysis_df = merge_by_common_ids(analysis_df, covariates)


#############################
# Cada PRS de ancestralidade entra como uma variável independente na
# regressão (da mesma forma que PC1, PC2, sexo, etc.), em vez de somar
# os 3 PRS em um único PRS_TOTAL.
#############################

prs_z_cols = []
for col in prs_cols:
    col_std = analysis_df[col].std()
    z_col = f"{col}_z"
    if pd.isna(col_std) or col_std == 0:
        analysis_df[z_col] = 0.0
        print(f"[AVISO] Desvio-padrão igual a zero para {col}; {z_col} definido como 0.")
    else:
        analysis_df[z_col] = (analysis_df[col] - analysis_df[col].mean()) / col_std
    prs_z_cols.append(z_col)

columns_to_use = prs_z_cols + covariates_to_include

n_antes = analysis_df.shape[0]
model_df = analysis_df.dropna(subset=columns_to_use + ["status"]).copy()
model_df["status"] = pd.to_numeric(model_df["status"], errors="coerce")
model_df = model_df.dropna(subset=["status"]).copy()
model_df = model_df.loc[model_df["status"].isin([0, 1])].copy()
n_removido = n_antes - model_df.shape[0]
if n_removido > 0:
    print(f"[AVISO] {n_removido} indivíduos removidos por NaN ou status inválido em PRS/status/covariáveis.")

y = model_df["status"].astype(int)

# ---- Modelo completo: PRS's + covariáveis ----
X_full = model_df[columns_to_use]
X_full = sm.add_constant(X_full)
modelo_full = sm.Logit(y, X_full).fit()

print(f"\n=== [PRS individuais + covariáveis] Resultado da regressão logística ===")
print(modelo_full.summary())

# ---- Modelo só com os PRS (sem covariáveis), usado para AUC e para o score combinado ----
X_prs_only = model_df[prs_z_cols]
X_prs_only = sm.add_constant(X_prs_only)
modelo_prs_only = sm.Logit(y, X_prs_only).fit()

print(f"\n=== [Somente os 3 PRS] Resultado da regressão logística ===")
print(modelo_prs_only.summary())

# pegando os valores principais de cada PRS separadamente, pra facilitar de olhar
prs_results = {}
for z_col in prs_z_cols:
    coef = modelo_full.params[z_col]
    pval = modelo_full.pvalues[z_col]
    ci_low, ci_high = modelo_full.conf_int().loc[z_col]

    odds_ratio = np.exp(coef)
    or_ci_low = np.exp(ci_low)
    or_ci_high = np.exp(ci_high)

    prs_results[z_col] = {
        "coef": coef,
        "pval": pval,
        "or": odds_ratio,
        "or_ci_low": or_ci_low,
        "or_ci_high": or_ci_high,
    }

    print(f"\n[{z_col}] (modelo com covariáveis)")
    print(f"Coeficiente (log-odds): {coef:.4f}")
    print(f"P-valor: {pval:.4g}")
    print(f"Odds Ratio (por DP de PRS): {odds_ratio:.3f} (IC95%: {or_ci_low:.3f} - {or_ci_high:.3f})")

y_true = y

# AUC só dos PRS's (poder discriminativo isolado dos 3 escores, sem covariáveis)
auc_prs = roc_auc_score(y_true, modelo_prs_only.predict(X_prs_only))

# AUC do modelo completo (PRS's + covariáveis) usando a probabilidade prevista
auc_full = roc_auc_score(y_true, modelo_full.predict(X_full))

print(f"\nAUC ({len(prs_z_cols)} PRS, sem covariáveis): {auc_prs:.4f}")
print(f"AUC ({len(prs_z_cols)} PRS + covariáveis): {auc_full:.4f}")

fpr_prs, tpr_prs, _ = roc_curve(y_true, modelo_prs_only.predict(X_prs_only))
fpr_full, tpr_full, _ = roc_curve(y_true, modelo_full.predict(X_full))

plt.figure(figsize=(6, 6))
plt.plot(fpr_prs, tpr_prs, label=f"{len(prs_z_cols)} PRS (AUC = {auc_prs:.3f})", color="darkorange")
if covariates_to_include:
    plt.plot(fpr_full, tpr_full, label=f"{len(prs_z_cols)} PRS + covariáveis (AUC = {auc_full:.3f})", color="firebrick")
plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Aleatório (AUC = 0.5)")
plt.xlabel("Taxa de Falso Positivo")
plt.ylabel("Taxa de Verdadeiro Positivo")
plt.title("Curva ROC — PRS (3 ancestralidades) predizendo status")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(
    f"{output_dir}/combined/{project_name}_PRS_total_roc_curve.png",
    dpi=150,
)
plt.close()  # libera a figura da memória antes da próxima ancestralidade no loop

# ---- distribuição do PRS: casos vs controles ----
# Como agora temos 3 PRS separados (e não mais um único PRS_TOTAL), usamos
# o preditor linear (log-odds) do modelo "só PRS" como um score combinado
# de 1 dimensão, apenas para fins de visualização da separação caso/controle.
# A regressão em si (acima) continua usando os 3 PRS como fatores independentes.

linear_pred = modelo_prs_only.predict(X_prs_only, linear=True)
lp_std = linear_pred.std()
if pd.isna(lp_std) or lp_std == 0:
    model_df["PRS_combined_z"] = 0.0
else:
    model_df["PRS_combined_z"] = (linear_pred - linear_pred.mean()) / lp_std

casos = model_df.loc[model_df["status"] == 1, "PRS_combined_z"]
controles = model_df.loc[model_df["status"] == 0, "PRS_combined_z"]

print(f"\nN casos: {len(casos)} | N controles: {len(controles)}")
print(f"Média score combinado casos: {casos.mean():.3f} | controles: {controles.mean():.3f}")

plt.figure(figsize=(7, 5))

x_range = np.linspace(model_df["PRS_combined_z"].min(), model_df["PRS_combined_z"].max(), 300)

kde_controles = gaussian_kde(controles)
kde_casos = gaussian_kde(casos)

plt.plot(x_range, kde_controles(x_range), color="steelblue", lw=2, label=f"Controle (n={len(controles)})")
plt.fill_between(x_range, kde_controles(x_range), alpha=0.2, color="steelblue")

plt.plot(x_range, kde_casos(x_range), color="firebrick", lw=2, label=f"Caso (n={len(casos)})")
plt.fill_between(x_range, kde_casos(x_range), alpha=0.2, color="firebrick")

plt.axvline(controles.mean(), color="steelblue", linestyle="--", lw=1.5,
            label=f"Média controle = {controles.mean():.3f}")
plt.axvline(casos.mean(), color="firebrick", linestyle="--", lw=1.5,
            label=f"Média caso = {casos.mean():.3f}")

plt.xlabel(f"Score combinado dos {len(prs_z_cols)} PRS (padronizado, z-score)")
plt.ylabel("Densidade")
plt.title(f"Distribuição do score combinado dos {len(prs_z_cols)} PRS: Caso vs Controle")
plt.legend()
plt.tight_layout()
plt.savefig(
    f"{output_dir}/combined/{project_name}_PRS_total_distribuicao_prs.png",
    dpi=150,
)
plt.close()

stat_u, p_u = mannwhitneyu(casos, controles, alternative="two-sided")
stat_t, p_t = ttest_ind(casos, controles, equal_var=False)

print(f"\nMann-Whitney U: p = {p_u:.4g}")
print(f"T-test (Welch): p = {p_t:.4g}")

run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
stats_output_path = os.path.join(
    output_dir,
    "combined",
    f"{project_name}_PRS_total_comparacao_estatistica_{run_timestamp}.txt"
)

summary_lines = [
    f"Ancestry: combined",
    f"Timestamp: {run_timestamp}",
    f"=== Resultado da regressão logística ({len(prs_z_cols)} PRS + covariáveis) ===",
    f"N casos: {len(casos)}",
    f"N controles: {len(controles)}",
    f"Média score combinado casos: {casos.mean():.4f}",
    f"Média score combinado controles: {controles.mean():.4f}",
    f"Mann-Whitney U p-value: {p_u:.4g}",
    f"T-test (Welch) p-value: {p_t:.4g}",
    "",
    "--- Coeficientes por PRS (modelo com covariáveis) ---",
]

for z_col, res in prs_results.items():
    summary_lines.extend([
        f"{z_col}:",
        f"  Coeficiente log-odds: {res['coef']:.4f}",
        f"  P-valor: {res['pval']:.4g}",
        f"  Odds Ratio: {res['or']:.3f}",
        f"  IC95% OR: {res['or_ci_low']:.3f} - {res['or_ci_high']:.3f}",
    ])

summary_lines.extend([
    "",
    f"AUC ({len(prs_z_cols)} PRS, sem covariáveis): {auc_prs:.4f}",
    f"AUC ({len(prs_z_cols)} PRS + covariáveis): {auc_full:.4f}",
])

write_summary_file(stats_output_path, summary_lines)

print(f"Resultados salvos em: {stats_output_path}")


print(covariates_to_include)

