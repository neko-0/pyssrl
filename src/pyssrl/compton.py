import numpy as np
import matplotlib.pyplot as plt

ELECTRON_REST_MASS = 5.11e5  # eV
ELECTRON_RADIUS_SQ = 0.07941  # barn
FINE_STRUCTURE_CONSTANT = 1.0 / 137.0


def klein_nishina(photon_e, angle):
    '''
    Klein-Nishina formula for differential cross section
    '''
    epsilon = photon_e / ELECTRON_REST_MASS
    inc_scatt_r = 1 / (1 + epsilon * (1 - np.cos(angle)))

    return (
        0.5
        * ELECTRON_RADIUS_SQ
        * inc_scatt_r**2
        * (inc_scatt_r + 1 / inc_scatt_r - np.sin(angle) ** 2)
    )


def recoil_electron_energy(photon_e, angle):
    '''
    photon_e (float) : incoming photon energy in eV.
    angle (float) : angle in radian.

    return:
        recoiled electron energy in eV.
    '''
    epsilon = photon_e / ELECTRON_REST_MASS
    l_ratio = 1.0 / (1.0 + epsilon * (1.0 - np.cos(angle)))
    scattered_e = photon_e * l_ratio
    return photon_e - scattered_e


def photoelectric_xsec(photon_e, Z=14):
    """
    This is simplified version
    """
    coeff = (
        16.0
        / 2.0
        * np.sqrt(2)
        * np.pi
        * ELECTRON_RADIUS_SQ
        * FINE_STRUCTURE_CONSTANT**4
    )
    k = photon_e / ELECTRON_REST_MASS
    return coeff * (Z**5) / (k**3.5)


def compton_xsec(photon_e, Z=14):
    """
    low energy Compton xsec (<100 keV)
    """
    k = photon_e / ELECTRON_REST_MASS
    coeff = Z * (8.0 / 3.0) * np.pi * ELECTRON_RADIUS_SQ
    a = 1.0 / (1.0 + 2 * k) ** 2
    return (
        coeff
        * a
        * (
            1
            + 2 * k
            + 6 / 5 * k**2
            - 0.5 * k**3
            + 2 / 7 * k**5
            + 8 / 105 * k**6
            + 4 / 105 * k**7
        )
    )


def pair_production_xsec(photon_e, Z=14):
    k = photon_e / ELECTRON_REST_MASS
    # the energy threshold should match the electron and positron rest mass
    if k < 2:
        return 0.0
    rho = (2 * k - 4) / (2 + k + 2 * np.sqrt(2 * k))
    coeff = Z**2 * FINE_STRUCTURE_CONSTANT * ELECTRON_RADIUS_SQ * (2 * np.pi / 3)
    # this is for low energies k < 4
    return (
        coeff
        * (1 - 2 / k) ** 3
        * (1 + 0.5 * rho + 23 / 40 * rho**2 + 11 / 60 * rho**3 + 29 / 960 * rho**4)
    )


def from_photoelectric(photon_e, nevents):
    compton = compton_xsec(photon_e)
    ph = photoelectric_xsec(photon_e)
    pair = pair_production_xsec(photon_e)
    # print(ph, compton, pair)

    # tot_xsec = ph + compton + pair
    tot_xsec = ph + compton

    ph_frac = ph / tot_xsec
    compton_frac = compton / tot_xsec
    pair_frac = pair / tot_xsec

    print(f"number of p-e events: {nevents}")
    print(f"photoelectric: {ph_frac} --> {int(ph_frac*nevents)}")
    print(f"compton: {compton_frac} --> {int(compton_frac*nevents)}")
    print(f"pair production: {pair_frac} --> {int(pair_frac*nevents)}")

    tot_evts = nevents / ph_frac
    compton_evts = tot_evts * compton_frac
    pair_evts = tot_evts * pair_frac

    return tot_evts, compton_evts


def _kn_recoil_sampling(photon_e, acceptance_angles=None):
    kappa = photon_e / ELECTRON_REST_MASS
    epsilon_0 = 1.0 / (1.0 + 2.0 * kappa)
    alpha1 = np.log(1 / epsilon_0)
    alpha2 = 0.5 * (1 - epsilon_0**2)
    alpha_sum = alpha1 + alpha2
    alpha1_frac = alpha1 / alpha_sum
    # alpha2_frac = alpha2 / alpha_sum
    rv_a = np.random.default_rng().uniform(0, 1)
    rv_b = np.random.default_rng().uniform(0, 1)
    rv_c = np.random.default_rng().uniform(0, 1)
    sign = np.random.choice([-1, 1])
    if rv_a < alpha1_frac:
        epsilon = np.exp(-np.log(1 / epsilon_0) * rv_b)
    else:
        epsilon = np.sqrt(epsilon_0**2 + (1 - epsilon_0**2) * rv_b)
    t = (1 - epsilon) / (kappa * epsilon)
    sin2_theta = t * (2 - t)
    g = 1 - epsilon * sin2_theta / (1 + epsilon**2)

    # angle = sign * np.arcsin(np.sqrt(sin2_theta))
    angle = sign * np.arccos(1 - ((1 / epsilon) - 1) / kappa)

    # if acceptance_angles:
    #     low, high = acceptance_angles
    #     low_epi = 1.0 / (1 + kappa*(1-np.cos(low)))
    #     high_epi = 1.0 / (1 + kappa*(1-np.cos(high)))
    #     if low_epi > high_epi:
    #         low_epi, high_epi = high_epi, low_epi
    #     if epsilon < low_epi or epsilon > high_epi:
    #         # print(angle*180/np.pi, low*180/np.pi, high*180/np.pi)
    #         return -1, angle

    if acceptance_angles:
        low, high = acceptance_angles
        if angle < low or angle > high:
            return -1, angle

    if rv_c < g:
        return photon_e - epsilon * photon_e, angle
    else:
        return -1, angle


def kn_recoil_sampling(
    photon_e, N=int(1e4), acceptance_angles=None, include_angles=False
):
    recoils = np.empty(N)
    angles = np.empty(N)
    for i in range(N):
        r_e, r_a = _kn_recoil_sampling(photon_e, acceptance_angles)
        recoils[i] = r_e
        angles[i] = r_a
    mask = recoils != -1
    if include_angles:
        return recoils[mask], angles[mask]
    else:
        return recoils[mask]


if __name__ == "__main__":
    # inc_e = 35e3
    # recoils = []
    # for i in range(1000):
    #     recoil_e = kn_recoil_sampling(inc_e)
    #     if recoil_e is None:
    #         continue
    #     recoils.append(recoil_e)

    # angles = np.arange(-np.pi, np.pi, 0.1)
    # xsec = klein_nishina(inc_e, angles)
    # recoil = recoil_electron_energy(inc_e, angles)

    # plt.plot(angles, xsec)
    # plt.plot(recoil, xsec)
    # plt.plot(angles*180/np.pi, recoil)

    # energies, angles = kn_recoil_sampling(35e3*2, include_angles=True,acceptance_angles=(1.021-0.066, 1.021+0.066))
    # plt.plot(angles*180/np.pi, energies)

    # angles = np.arange(0, 5, 0.1)
    # xsec = klein_nishina(60e3, angles)
    # recoil = recoil_electron_energy(60e3, angles)
    # plt.plot(recoil, xsec)

    # inc_e = 30e4
    # angles = np.arange(0, 5, 0.1)
    # xsec = diff_xsec(inc_e, angles)
    # recoil = recoil_e(inc_e, angles)
    # plt.plot(recoil/inc_e, xsec)

    # inc_e = 30e5
    # angles = np.arange(0, 5, 0.1)
    # xsec = diff_xsec(inc_e, angles)
    # recoil = recoil_e(inc_e, angles)
    # plt.plot(recoil/inc_e, xsec)
    # plt.ylim([0,0.1])

    # plt.show()
    def draw_energy_angles():
        inc_e_list = [70e3, 35e3, 30e3, 25e3, 10e3]
        # solid_angles = None
        solid_angles = (1.021 - 0.066, 1.021 + 0.066)

        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

        for inc_e in inc_e_list:
            energies, angles = kn_recoil_sampling(
                inc_e, include_angles=True, acceptance_angles=solid_angles
            )
            ax.scatter(angles, energies / 1e3, label=f"Inc. E={int(inc_e/1e3)}keV")

        ax.set_rticks([1, 2, 3, 4, 5])
        ax.grid(True)
        plt.legend()
        plt.show()

    def draw_xsec_angles():
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        inc_e_list = [10e6, 70e3, 35e3, 30e3, 25e3, 10e3]
        angles = np.arange(0, 2 * np.pi, 0.01)
        for inc_e in inc_e_list:
            xsec = klein_nishina(inc_e, angles)
            ax.scatter(angles, xsec, label=f"Inc. E={int(inc_e/1e3)}keV")
        # ax.set_rticks([1, 2, 3, 4, 5])
        ax.grid(True)
        plt.legend()
        plt.show()

    def draw_recoil_energy_xsec(photon_e):
        fig, ax = plt.subplots()
        angles = np.arange(0, 2 * np.pi, 0.01)
        for pe in photon_e:
            xsec = klein_nishina(pe, angles)
            recoil_e = recoil_electron_energy(pe, angles)
            ax.plot(recoil_e / 1e3, xsec, label=f"{int(pe/1e3)} keV")
        plt.legend()
        plt.show()

    # draw_energy_angles()
    # draw_xsec_angles()

    from_photoelectric(23e3, 10000)
    draw_recoil_energy_xsec([10e3, 20e3, 30e3, 60e3])
