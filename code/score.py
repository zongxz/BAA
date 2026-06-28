import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr



data = pd.read_csv("pred_age_meta.csv")

X = data.iloc[:, 16:].values.astype(float)
y = data["age_acceleration"].values



X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    random_state=42
)


scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)



X_train = np.expand_dims(X_train, axis=2)
X_test = np.expand_dims(X_test, axis=2)



model = Sequential()

# ---- Block 1 ----
model.add(Conv1D(32, 2, strides=1, activation='selu',
                 input_shape=(X_train.shape[1], 1)))
model.add(Conv1D(32, 2, strides=1, activation='selu'))
model.add(MaxPooling1D(2))

# ---- Block 2 ----
model.add(Conv1D(32, 2, strides=1, activation='selu'))
model.add(Conv1D(32, 2, strides=1, activation='selu'))
model.add(MaxPooling1D(2))

# ---- Block 3 ----
model.add(Conv1D(32, 2, strides=1, activation='selu'))
model.add(Conv1D(32, 2, strides=1, activation='selu'))
model.add(MaxPooling1D(2))

# ---- Block 4 (last conv layers to reach 8 total) ----
model.add(Conv1D(32, 2, strides=1, activation='selu'))
model.add(Conv1D(32, 2, strides=1, activation='selu'))

# ---- Head ----
model.add(Flatten())
model.add(Dense(32, activation='selu'))
model.add(Dense(32, activation='selu'))
model.add(Dense(1, activation='linear'))


# Compile


model.compile(
    optimizer=keras.optimizers.Adam(learning_rate=1e-4),
    loss='mse',
    metrics=['mae']
)




callback = keras.callbacks.EarlyStopping(
    monitor='val_loss',
    patience=20,
    restore_best_weights=True
)



history = model.fit(
    X_train,
    y_train,
    validation_split=0.2,
    epochs=500,
    batch_size=1024,
    callbacks=[callback],
    verbose=1
)



y_pred = model.predict(X_test).flatten()



mae = mean_absolute_error(y_test, y_pred)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
r2 = r2_score(y_test, y_pred)
pearson_r = pearsonr(y_test, y_pred)[0]

print("MAE:", mae)
print("RMSE:", rmse)
print("R2:", r2)
print("Pearson r:", pearson_r)