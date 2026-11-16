#ifndef CECE_MEGAN_HISTORY_HPP
#define CECE_MEGAN_HISTORY_HPP

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <stdexcept>
#include <string>

#include "cece/cece_logger.hpp"

namespace cece {

using MeganHistoryInputView = Kokkos::View<const double***, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace>;
using MeganHistoryView = Kokkos::View<double***, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace>;

/**
 * @brief Per-cell running histories used by the HEMCO MEGAN implementation.
 *
 * State is cold-started on the first run and updated after emissions with the
 * HEMCO exponential running-mean weights (dt / tau). Restart I/O is deliberately
 * left to a future integration layer.
 *
 * temperature_12hour mirrors HEMCO's T_LAST24H/T_PREVDAY. HEMCO 3.12.1 passes it
 * to GET_GAMMA_T_LD but documents that argument as unused, so no emission term
 * reads it; it is kept for diagnostics and comparison.
 */
class MeganHistoryState {
   public:
    MeganHistoryView temperature_5day;
    MeganHistoryView temperature_12hour;
    MeganHistoryView par_direct_5day;
    MeganHistoryView par_diffuse_5day;
    MeganHistoryView lai_12hour;

    void InitializeIfNeeded(int nx, int ny, const MeganHistoryInputView& lai, double temperature_5day_init, double temperature_12hour_init,
                            double par_direct_init, double par_diffuse_init) {
        if (initialized_) {
            if (nx != nx_ || ny != ny_) {
                throw std::runtime_error("MeganHistoryState: grid shape changed after history initialization");
            }
            return;
        }

        nx_ = nx;
        ny_ = ny;
        temperature_5day = MeganHistoryView("megan_temperature_5day", nx, ny, 1);
        temperature_12hour = MeganHistoryView("megan_temperature_12hour", nx, ny, 1);
        par_direct_5day = MeganHistoryView("megan_par_direct_5day", nx, ny, 1);
        par_diffuse_5day = MeganHistoryView("megan_par_diffuse_5day", nx, ny, 1);
        lai_12hour = MeganHistoryView("megan_lai_12hour", nx, ny, 1);

        auto t5 = temperature_5day;
        auto t12 = temperature_12hour;
        auto pd = par_direct_5day;
        auto pf = par_diffuse_5day;
        auto l12 = lai_12hour;
        Kokkos::parallel_for(
            "MeganHistoryState_ColdStart", Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nx, ny}),
            KOKKOS_LAMBDA(int i, int j) {
                t5(i, j, 0) = temperature_5day_init;
                t12(i, j, 0) = temperature_12hour_init;
                pd(i, j, 0) = par_direct_init;
                pf(i, j, 0) = par_diffuse_init;
                l12(i, j, 0) = lai(i, j, 0);
            });
        Kokkos::fence();
        initialized_ = true;
    }

    /// Evolving histories need the scheduler's step length and calendar day.
    static void RequireClock(const std::string& who, int dt_seconds, int day_of_year) {
        if (dt_seconds <= 0 || day_of_year <= 0) {
            throw std::runtime_error(who +
                                     ": evolving history requires a scheduler run context with positive dt_seconds and day_of_year; "
                                     "run through the CECE clock or set history_mode: fixed");
        }
    }

    /// Update the running means with weight dt / tau, clamped to 1 (warning once) when dt exceeds tau.
    void Update(const MeganHistoryInputView& temperature, const MeganHistoryInputView& par_direct, const MeganHistoryInputView& par_diffuse,
                const MeganHistoryInputView& lai, int dt_seconds) {
        if (!initialized_) {
            throw std::runtime_error("MeganHistoryState: cannot update before cold-start initialization");
        }
        if (dt_seconds <= 0) {
            throw std::runtime_error("MeganHistoryState: dt_seconds must be positive");
        }
        if (dt_seconds > kTauShortSeconds && !clamp_warned_) {
            CECE_LOG_WARNING("MeganHistoryState: dt of " + std::to_string(dt_seconds) +
                             " s exceeds the 12-hour history e-folding time; running-mean weights are clamped to 1 "
                             "(history equals the current sample). Further occurrences are not reported.");
            clamp_warned_ = true;
        }

        const double daily_new_fraction = std::min(1.0, static_cast<double>(dt_seconds) / kTauTemperatureParSeconds);
        const double daily_old_fraction = 1.0 - daily_new_fraction;
        const double short_new_fraction = std::min(1.0, static_cast<double>(dt_seconds) / kTauShortSeconds);
        const double short_old_fraction = 1.0 - short_new_fraction;

        auto t5 = temperature_5day;
        auto t12 = temperature_12hour;
        auto pd = par_direct_5day;
        auto pf = par_diffuse_5day;
        auto l12 = lai_12hour;
        Kokkos::parallel_for(
            "MeganHistoryState_Update", Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nx_, ny_}),
            KOKKOS_LAMBDA(int i, int j) {
                t5(i, j, 0) = daily_old_fraction * t5(i, j, 0) + daily_new_fraction * temperature(i, j, 0);
                pd(i, j, 0) = daily_old_fraction * pd(i, j, 0) + daily_new_fraction * par_direct(i, j, 0);
                pf(i, j, 0) = daily_old_fraction * pf(i, j, 0) + daily_new_fraction * par_diffuse(i, j, 0);
                t12(i, j, 0) = short_old_fraction * t12(i, j, 0) + short_new_fraction * temperature(i, j, 0);
                l12(i, j, 0) = short_old_fraction * l12(i, j, 0) + short_new_fraction * lai(i, j, 0);
            });
        Kokkos::fence();
    }

   private:
    static constexpr int kTauTemperatureParSeconds = 5 * 24 * 60 * 60;
    static constexpr int kTauShortSeconds = 12 * 60 * 60;

    int nx_ = 0;
    int ny_ = 0;
    bool initialized_ = false;
    bool clamp_warned_ = false;
};

}  // namespace cece

#endif  // CECE_MEGAN_HISTORY_HPP
