import numpy as np
import pandas as pd
from mcarma import ObservationData, fit

# Read the same observations used by Rdrw
df = pd.read_csv("mGP_ex1_fiveband.csv")

print("Observations:", len(df))
print("Counts by band:")
print(df.groupby("band").size())

# Prepare input arrays
t = df["time"].to_numpy(dtype=float)
y = df["value"].to_numpy(dtype=float)
err = df["error"].to_numpy(dtype=float)
band = df["band"].to_numpy(dtype=int)

# Construct observation data
data = ObservationData(
    t_obs=t,
    y_obs=y,
    band=band,
    R=err**2,
    d=5
)

# Fit MCARMA(1,0)
res = fit(
    data,
    p=1,
    q=0,
    n_restarts=5,
    use_jax_grad=True
)

print("\nMCARMA(1,0) results")
print("--------------------")
print("Unpenalized logLik:", res["loglik_pure"])
print("Penalized logLik:", res["loglik"])
print("Estimated parameters:")
print(res["theta"])

