# systems_biology_models.R
# Author: Fares Ibrahim
# Description: Collection of systems biology models in R (ODEs, curve fitting, kinetics)

# Required Libraries
library(deSolve)      # For solving ODEs
library(ggplot2)      # For visualization
library(nlstools)     # For nonlinear fitting
library(stats)        # Built-in fitting functions

#----------------------------#
# 1. GUT–BLOOD ODE MODEL
#----------------------------#

# Parameters
half_life <- 5        # Vancomycin half-life (hours)
k <- log(2)/half_life # Elimination rate constant
a <- 0.9              # Absorption rate (/h)
G <- 1000             # Initial gut concentration
B <- 0                # Initial blood concentration

# Initial Conditions
init <- c(G = G, B = B)
t <- seq(0, 40, 0.01)

# ODE Function
gut_blood_ode <- function(t, state, parms) {
  with(as.list(state), {
    dG <- -a * G
    dB <- a * G - k * B
    list(c(dG, dB))
  })
}

# Solve
out <- ode(y = init, times = t, func = gut_blood_ode, parms = NULL)
out <- as.data.frame(out)

# Plot
plot(out$time, out$G, type = "l", col = "blue", ylab = "Concentration", xlab = "Time (h)")
lines(out$time, out$B, col = "red")
legend("topright", legend = c("Gut", "Blood"), col = c("blue", "red"), lty = 1)


#----------------------------#
# 2. ORAL PK MODEL (Analytical)
#----------------------------#

ka <- 0.9
dose <- 1000
ke <- log(2)/half_life
time <- seq(0, 24, 0.01)

Ct <- (dose * ka / (ka - ke)) * (exp(-ke * time) - exp(-ka * time))

plot(time, Ct, type = "l", col = "darkgreen", xlab = "Time (h)", ylab = "Plasma Concentration")

t_max <- time[which.max(Ct)]
cat("Time of peak concentration (Tmax):", round(t_max, 2), "h\n")


#----------------------------#
# 3. LIGAND–RECEPTOR BINDING
#----------------------------#

binding_model <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dC <- kon * (Rtot - C) * L - koff * C
    list(c(dC))
  })
}

params <- c(kon = 0.6, koff = 0.01, Rtot = 2000, L = 0.0167)
state <- c(C = 0)
times <- seq(0, 300, 0.1)

binding_out <- ode(y = state, times = times, func = binding_model, parms = params)
binding_out <- as.data.frame(binding_out)

# Plot
plot(binding_out$time, binding_out$C, type = "l", col = "blue", xlab = "Time (s)", ylab = "Bound Receptors")

# Receptors bound at 3s
receptors_at_3 <- binding_out$C[which.min(abs(binding_out$time - 3))]
abline(h = receptors_at_3, col = "blue", lty = 2)
abline(v = 3, col = "red", lwd = 2)
text(0, receptors_at_3, labels = round(receptors_at_3, 2), pos = 4)

# Steady-state
final <- tail(binding_out$C, 1)
threshold <- 0.95 * final
ss_time <- min(binding_out$time[binding_out$C >= threshold])
abline(v = ss_time, col = "green", lty = 2)
text(ss_time, 0.9 * max(binding_out$C), labels = paste("Steady at", round(ss_time, 2)), col = "green")


#----------------------------#
# 4. ENZYME KINETICS
#----------------------------#

NADH <- c(0, 0.075, 0.150, 0.225, 0.300, 0.375)
Abs <- c(0, 0.5065, 1.1743, 1.2, 1.48, 1.7)

df1 <- data.frame(NADH, Abs)
fit <- lm(Abs ~ NADH - 1, data = df1)
slope <- coef(fit)

plot(df1$NADH, df1$Abs, xlab = "NADH (mM)", ylab = "Absorbance")
abline(fit, col = "red")

# Timepoint data
pyruvate <- c(0.15, 0.3, 0.6, 1.2, 1.8)
Abs_t10 <- c(1.129, 0.991, 0.972, 0.746, 0.879)
Abs_t20 <- c(1.087, 0.915, 0.872, 0.628, 0.760)

NADH_t10 <- Abs_t10 / slope
NADH_t20 <- Abs_t20 / slope

V <- (NADH_t10 - NADH_t20) / 10
df2 <- data.frame(pyruvate, V)

# Fit MM model
mm_fit <- nls(V ~ Vmax * pyruvate / (Km + pyruvate), data = df2,
              start = list(Vmax = 0.008, Km = 0.5))

summary(mm_fit)

# CI Calculation
coefs <- summary(mm_fit)$coefficients[, 1]
SEs <- summary(mm_fit)$coefficients[, 2]
CI <- cbind(coefs - qt(0.975, df = 3) * SEs,
            coefs + qt(0.975, df = 3) * SEs)
colnames(CI) <- c("Lower", "Upper")
CI

# Plot
plot(pyruvate, V, ylim = c(0, 0.0035), xlim = c(0, 2),
     xlab = "Pyruvate (mM)", ylab = "Velocity (mM/s)")
curve(coefs[1] * x / (coefs[2] + x), col = "red", add = TRUE)


#----------------------------#
# 5. SIRS MODEL
#----------------------------#

sirs <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dS <- -beta * I * S / N + delta * R
    dI <- beta * I * S / N - gamma * I
    dR <- gamma * I - delta * R
    list(c(dS, dI, dR))
  })
}

params <- c(beta = 0.4, gamma = 0.2, delta = 0.005)
state <- c(S = 995, I = 5, R = 0)
N <- 1000
times <- seq(0, 1000, 1)

sirs_out <- ode(y = state, times = times, func = sirs, parms = params)
sirs_out <- as.data.frame(sirs_out)

matplot(sirs_out$time, sirs_out[, c("S", "I", "R")], type = "l", lty = 1,
        col = c("blue", "red", "green"), ylab = "Individuals", xlab = "Time")
legend("right", legend = c("Susceptible", "Infectious", "Recovered"),
       col = c("blue", "red", "green"), lty = 1)
