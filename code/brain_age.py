import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

import statsmodels.api as sm


data = pd.read_csv("brain_mri_i2.csv")

y = data["age_i2"]
X = data.drop(columns=["age_i2"])


X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    random_state=42
)


model = Pipeline([
    ("scaler", StandardScaler()),
    ("lasso", LassoCV(cv=5, random_state=42, max_iter=10000))
])

model.fit(X_train, y_train)


pred_train = model.predict(X_train)
pred_test = model.predict(X_test)



def evaluate(y_true, y_pred, name="set"):
    r2 = r2_score(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))

    print(f"\n{name}")
    print(f"R2: {r2:.4f}")
    print(f"MAE: {mae:.4f}")
    print(f"RMSE: {rmse:.4f}")

evaluate(y_train, pred_train, "Train")
evaluate(y_test, pred_test, "Test")


pred_result = pd.DataFrame({
    "age": y_test.values,
    "pred_age": pred_test
})

pred_result["brain_age_gap"] = pred_result["pred_age"] - pred_result["age"]


bias_model = sm.OLS(
    pred_result["pred_age"],
    sm.add_constant(pred_result["age"])
).fit()

alpha, beta = bias_model.params

pred_result["cal_pred_age"] = pred_result["pred_age"] + pred_result["age"] - (alpha + beta * pred_result["age"])


pred_result["cal_brain_age_gap"] = pred_result["cal_pred_age"] - pred_result["age"]


evaluate(pred_result["age"], pred_result["pred_age"], "Before correction")
evaluate(pred_result["age"], pred_result["cal_pred_age"], "After correction")



lasso = model.named_steps["lasso"]

importance = pd.DataFrame({
    "feature": X.columns,
    "coef": lasso.coef_
})

importance["abs_coef"] = np.abs(importance["coef"])
importance = importance.sort_values("abs_coef", ascending=False)

importance.to_csv("feature_importance.csv", index=False)
