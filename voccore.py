from dateutil import rrule, parser
import numpy as np
from datetime import datetime             # for manipulating dates
import matplotlib.pyplot as plt           # for plotting
from scipy.stats import norm, uniform, lognorm, triang
from scipy.ndimage import gaussian_filter1d
import scipy.special as sps
from tqdm import tqdm
import time
from datetime import  timedelta

from scipy.integrate import odeint

def sample_parameter_space(config):
    """
    Sample the parameter space defined in the config file. Take n_samples
    :param config: Dict containing configuration options
    :return: Sampled parameter space (dict of iterables)
    """
    n_samples = config['n_samples']

    date1 = '2021-01-07'
    date2 = '2021-06-01'
    datesx = list(rrule.rrule(rrule.DAILY, dtstart=parser.parse(date1), until=parser.parse(date2)))

    x = np.array(datesx)
    n_days = np.size(datesx)

    sample_space = dict()
    variable_param_names = ['r_lockdownscale', 'k_voc', 'f0_voc','vac_ratio',
                           'startinfperc','startsusceptible', 'ts','ve_trans',
                            've_immune_o', 've_immune_hosp_o',
                            've_vac_d', 've_vac_o',
                            've_booster_d', 've_booster_o',
                            've_vac_hosp_d', 've_vac_hosp_o',
                            've_booster_hosp_d', 've_booster_hosp_o',
                            're0_voc', 'seir_tlatent','seir_tinf','tshift']
    locked_param_names = ['calcmethod', 'startdate', 'enddate',
                          'p_vac_start', 'p_booster_start',
                          'r_lockdowndayx', 'r_lockdown',
                          'boosterdayx', 'boosterdayn',
                          'demage', 'demn',
                          'agegroups_p_vac', 'w_hosp_age']
    np.random.seed(1)
    varying_parameters = []
    for param_name in variable_param_names:
        try:
            config_name = config[param_name][0].lower()
            varying_parameters.append(param_name)
            if config_name == 'normal':
                mean = config[param_name][1]
                std_dev = config[param_name][2]
                sample_space[param_name] = norm(loc=mean, scale=std_dev).rvs(size=n_samples)
            elif config_name == 'uniform':
                min = config[param_name][1]
                max = config[param_name][2]
                sample_space[param_name] = uniform(loc=min, scale=max - min).rvs(size=n_samples)
            elif config_name == 'lognormal':
                exp_mu = np.exp(config[param_name][1])
                sigma = config[param_name][2]
                sample_space[param_name] = lognorm(s=sigma, scale=exp_mu).rvs(size=n_samples)
            elif config_name == 'triangular':
                min = config[param_name][1]
                peak = config[param_name][2]
                max = config[param_name][3]
                sample_space[param_name] = triang(loc=min, scale=max - min, c=(peak - min) / (max - min)).rvs(
                    size=n_samples)
        except TypeError:
            sample_space[param_name] = np.full(n_samples, fill_value=float(config[param_name]))

    for param_name in locked_param_names:
        sample_space[param_name] = [config[param_name] for _ in range(n_samples)]

    sample_space['time'] = np.linspace(0, n_days, n_days + 1).astype(float)


    return sample_space, varying_parameters

def run_models(config):
    n_samples = config['n_samples']
    p, varying_param = sample_parameter_space(config)
    print('Running {} models'.format(n_samples))
    print('Varying: {}'.format(', '.join(varying_param)))
    time.sleep(0.1)  # For pretty printing

    results = []
    for a in tqdm(range(n_samples)):
        out = single_run (p['r_lockdowndayx'][a], p['r_lockdown'][a], p['r_lockdownscale'][a],
                                 p['k_voc'][a], p['f0_voc'][a], p['vac_ratio'][a],
                                 p['startdate'][a], p['enddate'][a],
                                 p['boosterdayx'][a], p['boosterdayn'][a],
                                 p['p_vac_start'][a],
                                 p['p_booster_start'][a],
                                 p['demage'][a], p['demn'][a],
                                 p['agegroups_p_vac'][a], p['w_hosp_age'][a],
                                 p['startinfperc'][a], p['startsusceptible'][a], p['ts'][a],
                                 ve_trans=p['ve_trans'][a],
                                 ve_immune_o=p['ve_immune_o'][a],
                                 ve_immune_hosp_o=p['ve_immune_hosp_o'][a],
                                 ve_vac_d=p['ve_vac_d'][a],
                                 ve_vac_o=p['ve_vac_o'][a],
                                 ve_booster_d=p['ve_booster_d'][a],
                                 ve_booster_o=p['ve_booster_o'][a],
                                 ve_vac_hosp_d=p['ve_vac_hosp_d'][a],
                                 ve_vac_hosp_o=p['ve_vac_hosp_o'][a],
                                 ve_booster_hosp_d=p['ve_booster_hosp_d'][a],
                                 ve_booster_hosp_o=p['ve_booster_hosp_o'][a],
                                 re0_voc=p['re0_voc'][a],
                                 calcmethod=p['calcmethod'][a],
                                 seir_tlatent=p['seir_tlatent'][a],
                                 seir_tinf=p['seir_tinf'][a],
                                 tshift=p['tshift'][a]
                             )

        results.append(out)

    return p, results


def get_expectedk_segregated (pu, pv, pb, S, ve_trans, ve_vac_d, ve_booster_d, ve_vac_o, ve_booster_o, ve_immune_o, R0o_R0d, ts):
    Ud = (1-S)*(1-ve_immune_o)/S

    Ct_o =  1
    Ct_d =  1

    Id = 1
    Io = (pu +Ud *(1 - ve_immune_o)*(1-ve_trans) +pv * (1 - ve_vac_o)*(1-ve_trans)  + pb * (1 - ve_booster_o)*(1-ve_trans) ) / (
                pu + pv * (1 - ve_vac_d)*(1-ve_trans)  + pb * (1 - ve_booster_d)*(1-ve_trans) )

    Rtd = Ct_d * Id

    gamma = 1 / ts
    #beta = np.exp(k) - 1 + gamma
    #rtratioref = beta / gamma
    # R0o_R0d = rtratioref * Ct_d * Id / (Ct_o * Io)
    rtratioref  = R0o_R0d * (Ct_o * Io)/ (Ct_d * Id)
    beta =   rtratioref * gamma
    expk = beta +1 - gamma
    k = np.log(expk)

    return k

def get_expectedk_mixed (pu, pv, pb, S, ve_trans, ve_vac_d, ve_booster_d, ve_vac_o, ve_booster_o, ve_immune_o, R0o_R0d, ts):
    Ud = (1-S)*(1-ve_immune_o)/S
    Cv =  (pu + (pv + pb) * (1 - ve_trans))
    Ct_o =  ((pu + (pv + pb + Ud) * (1 - ve_trans)) / (pu + (pv + pb + Ud))) / Cv
    Ct_d =  (pu + (pv + pb) * (1 - ve_trans)) / Cv

    Id = 1
    Io = (pu + Ud*(1 - ve_immune_o) +  pv  * (1 - ve_vac_o) + pb * (1 - ve_booster_o)) / (
                pu + pv * (1 - ve_vac_d) + pb * (1 - ve_booster_d))

    Rtd = Ct_d * Id

    gamma = 1 / ts
    #beta = np.exp(k) - 1 + gamma
    #rtratioref = beta / gamma
    # R0o_R0d = rtratioref * Ct_d * Id / (Ct_o * Io)
    rtratioref  = R0o_R0d * (Ct_o * Io)/ (Ct_d * Id)
    beta =   rtratioref * gamma
    expk = beta +1 - gamma
    k = np.log(expk)
    return k


def getHosp(inf, ft, demage, demn, agegroups_p_vac, w_hosp_age, p_booster_start, p_vac_start, dp_b,
            Ud, ve_immune_hosp_o, ve_vac_hosp_d, ve_vac_hosp_o, ve_booster_hosp_d, ve_booster_hosp_o):
    # in each agegroups assess the expected vaccination levels
    demtotal = demn[-1]
    agegroupsn = np.diff(demn)

    #establish which fractions inside the different age groups are vaccinated or boostered
    agegroups = demage[0:-1]
    # target vaccination rate in age groups (use this as maximum to fill with available vaccins and boosters,
    # start from oldest group
    agegroupspvac = agegroups_p_vac
    # the weight of the agegroup
    w = agegroupsn / demtotal
    #  the weight fraction vaccins relative to total population
    wpvac = w * agegroupspvac
    agegroup_booster = np.empty([np.size(agegroups), np.size(inf)])
    agegroup_vac = np.empty([np.size(agegroups), np.size(inf)])

    # fill the vaccins and boosters in the different age groups
    for j, index in enumerate(inf):
        p_booster = p_booster_start + dp_b[j]
        p_vac =  p_vac_start - dp_b[j]
        for i, agemax in enumerate(agegroupsn):
            agegroup_booster[i][j] = min(wpvac[i], p_booster)
            p_booster -= agegroup_booster[i][j]
            agegroup_vac[i][j] = min((wpvac[i]-agegroup_booster[i][j]), p_vac)
            p_vac -= agegroup_vac[i][j]


    # Now determine the hospitalization pressure from infected
    # first establish the reference weights of each age group on the hospitalization, which equals w*CHRd0
    w_hosp = w_hosp_age/np.sum(w_hosp_age)

    hosp = inf*0
    for i, agemax in enumerate(agegroupsn):
        pv_iage = agegroup_vac[i]/w[i]
        pb_iage = agegroup_booster[i]/w[i]
        pu_iage = 1 - pv_iage -pb_iage

        hdelta = (pu_iage + pv_iage * (1 - ve_vac_hosp_d) + pb_iage * (1 - ve_booster_hosp_d))
        homicron =  ((1-ve_immune_hosp_o) *pu_iage + (pv_iage +Ud) * (1 - ve_vac_hosp_o) + pb_iage * (1 - ve_booster_hosp_o))/ (pu_iage + pv_iage + pb_iage + Ud)

        hosp_agegroup = (1-ft) *hdelta + ft * homicron
        hosp += hosp_agegroup *w_hosp[i]

    hosp *= inf
    hosp = hosp/hosp[0]
    return hosp

def single_run (r_lockdowndayx, r_lockdownval, r_lockdownscale,
               k_voc,  f0_voc, vac_ratio, startdate, enddate, boosterdayx, boosterdayn,
               p_vac_start,  p_booster_start,  demage, demn,
               agegroups_p_vac, w_hosp_age,
               startinfperc, startsusceptible,  ts,
               ve_trans = 0.5,
               ve_immune_o = 0.2, ve_immune_hosp_o = 0.5,
               ve_vac_d = 0.6, ve_vac_o=0.2,
               ve_booster_d=0.95, ve_booster_o=0.75,
               ve_vac_hosp_d= 0.95, ve_vac_hosp_o=0.48,
               ve_booster_hosp_d=0.975, ve_booster_hosp_o=0.9,
               re0_voc=-1,
               calcmethod=1,
               seir_tlatent=2,
               seir_tinf=5, tshift=14):
    """
         construct Grid instance from fname and format

     :param  r_lockdowndayx: days after start specifying r_lockdown (for delta)
     :param  r_lockdownval: the R(0) value at r_lockdowndayx  (for delta)
     :param: r(0) scaling value
     :param  k_voc : the daily growth rate (natural logarithm of the VOC (omicron))
     :param  f0_voc  log10 value of f(t) at starting time
     :param  vac_ratio: the  booster speed as fraction of boosterdayx, boosterdayn
     :param  startdate  starting date
     :param  enddate  end date

      :param  boosterdayx: array with days from start
     :param  boosterdayn: number of boosters planned at boosterdayx
     :param  p_vac_start:  fraction of population vaccinated at start of simulation
     :param  p_booster_start:  fraction of population boostered at start of simulation

     :param  demage :  array with age cohort limits in the demography (descending, starting from ca 90 to 0)
     :param  demn : number of people in the age cohort
     :param  startinfperc: starting infection fraction of population (daily infected * ts)
     :param  startsusceptable: starting susceptable fraction of population (taking into account passed infections)
     :param  ts : serial interval time
     :param  ve_trans  vaccin and booster efficy against transmission
     :param  ve_immune_o  fraction in Recovered protected from omicron
     :param  ve_immune_hosp_o  unvaccinated protected againtst hospitalization (relative to delta)
     :param  ve_vac_d:  Vaccin efficacy against infections for delta
     :param  ve_vac_o   Vaccin efficacy against infections for omicron
     :param  ve_booster_d:  booster efficacy against infections for delta
     :param  ve_booster_o   booster efficacy against infections for omicron
     :param  ve_vac_hosp_d:  Vaccin  efficacy against hospitalization for delta
     :param  ve_vac_hosp_o   Vaccin (and delta infected loosing immunity) efficacy against hospitalization for omicron
     :param  ve_booster_hosp_d:  booster efficacy against hospitalization for delta
     :param  ve_booster_hosp_o   booster efficacy against hospitalization for omicron
     :param  re0_voc:  the ratio of the basic reproduction number of the VOC omicron  to prevailing (delta), if -1 it will estimate it from the vaccination state and k value
     :param  calcmethod:  calculation method 1: SIR, 2 SEIR, for SEIR seir_tlatent is latency time,  ts is infectious time
     :param  seit_tlatent: latency time (in days) for SEIR method (if used)
     :param  seit_tinf: infecteous time (in days_ for SEIR method (if used) (in order to arrrive at 4.8 days ts: tlatent and tinf van be set to 2 and 5 days)
     :param  tshift:  shift of reference Rt values reference dates as they start the be felt in registered cases +tshift days


     """

    iextra_conv = 20
    date1 = startdate
    date2 = enddate
    datesx = list(rrule.rrule(rrule.DAILY, dtstart=parser.parse(date1), until=parser.parse(date2)))

    date2ext = datesx[-1]+ timedelta(days=iextra_conv)
    datesx = list(rrule.rrule(rrule.DAILY, dtstart=parser.parse(date1), until=date2ext))

    x = np.array(datesx)
    indx = np.arange(x.size)

    ft = np.empty([indx.size])


    k = k_voc







    p = 10**f0_voc

    ft = p * np.exp(k* (indx)) / (p * np.exp(k * (indx)) + (1 - p))


    demtotal = demn[-1]

    # number of boosters
    dembooster = np.interp(indx-tshift, boosterdayx, boosterdayn)
    dembooster *=vac_ratio

    # you cannot booster more than vaccinated
    dembooster = np.minimum(dembooster, demtotal*(p_vac_start))

    dp_b = dembooster/demtotal

    # unvaccinated kept fixed
    pu = ft*0 + (1.0 -  p_vac_start - p_booster_start)
    # create time series of pv and pb
    pv = ft * 0 + p_vac_start - dp_b
    pb = ft * 0 + p_booster_start + dp_b



    Ud = (1-startsusceptible)*(1-ve_immune_o) / startsusceptible




    # eq 2 correction factors
    Cv =  (pu[0] + (p_vac_start + p_booster_start) * (1 - ve_trans))
    Ct_o = ft * 0 + ((pu + (pv + pb + Ud) * (1 - ve_trans)) / (pu + (pv + pb + Ud))) / Cv
    Ct_d = ft * 0 + (pu + (pv + pb) * (1 - ve_trans)) / Cv

    Id = ft * 0 + (pu + pv * (1 - ve_vac_d) + pb * (1 - ve_booster_d)) / (pu[0] + p_vac_start * (1 - ve_vac_d) + p_booster_start * (1 - ve_booster_d))
    Io = ft * 0 + (pu + Ud *(1-ve_immune_o) + pv* (1 - ve_vac_o) + pb * (1 - ve_booster_o)) / (pu[0] + (p_vac_start) * (1 - ve_vac_d) + p_booster_start * (1 - ve_booster_d))


    Rtd = Ct_d * Id

    gamma = 1 / ts
    beta = np.exp(k) - 1 + gamma
    rtratioref = beta / gamma
    if (re0_voc<0):
        R0o_R0d = rtratioref * Ct_d[0]*Id[0]/(Ct_o[0]*Io[0])
    else:
        R0o_R0d = re0_voc


    Rto=   R0o_R0d* Ct_o * Io


    rtratio_voc = Rto/Rtd

    #rtratio_voc, R0o_R0d = get_rtratio_voc(pu, pv, pb, Ud, Vo, Bo, Vd, Bd, p_vac_start, p_booster_start, k, ts)





    Rt = np.empty([indx.size])



    Rt = (1 - ft) *  Rtd+ ft * Rto




    # now correct properly for R(0) and measures for additional lockdown
    r_lockdown = np.interp(indx-tshift, r_lockdowndayx, r_lockdownval)
    Rt *= r_lockdown/Rt[0]



    inf = np.empty([indx.size])
    infected = startinfperc*1e-2
    #infected = 10e3 * ts / 17e6

    # should be corrected to get recovered back otherwise underestimate
    startsusceptible += Ud*startsusceptible
    sirinit = [startsusceptible - infected, infected, 0.0]




    if (calcmethod==1):
        inf = base_sir_model(sirinit, Rt, ts)
    else:
        inf = base_seir_model(sirinit, Rt, k,ts, seir_tlatent, seir_tinf)


    hosp = getHosp(inf, ft, demage, demn, agegroups_p_vac, w_hosp_age, p_booster_start, p_vac_start, dp_b,
            Ud, ve_immune_hosp_o,ve_vac_hosp_d, ve_vac_hosp_o, ve_booster_hosp_d, ve_booster_hosp_o)


    # next convolve the time series

    gamma_mean = 7 # 7
    gamma_stdev = 4 # 4

    inf = lognormal_smooth_shift_convolve(inf, gamma_mean, gamma_stdev, scale=1.0, fillval=inf[0])
    inf = lognormal_smooth_shift_convolve(inf, -gamma_mean, 0, scale=1.0, fillval=inf[0])

    gamma_mean_hosp = gamma_mean
    gamma_stdev_hosp =gamma_stdev

    hosp = lognormal_smooth_shift_convolve(hosp, gamma_mean_hosp, gamma_stdev_hosp, scale=1.0, fillval=inf[0])
    hosp = lognormal_smooth_shift_convolve(hosp, -gamma_mean, 0, scale=1.0, fillval=inf[0])






    #hosp = getHosp(inf, ft, demage, demn, agegroups_p_vac, w_hosp_age, p_booster_start, p_vac_start, dp_b,
    #       Ud, ve_immune_hosp_o,ve_vac_hosp_d, ve_vac_hosp_o, ve_booster_hosp_d, ve_booster_hosp_o)

    inf = inf / inf[0]

    Rt = getRt(inf, ts)

    infpeak = getpeaks(inf)





    hosppeak = getpeaks(hosp)

    # trim the results
    indx = indx[:-iextra_conv]
    inf = inf[:-iextra_conv]
    ft = ft[:-iextra_conv]
    Rt = Rt[:-iextra_conv]
    r_lockdown = r_lockdown[:-iextra_conv]
    rtratio_voc =rtratio_voc[:-iextra_conv]
    pu = pu[:-iextra_conv]
    pv = pv[:-iextra_conv]
    pb = pb[:-iextra_conv]
    hosp = hosp[:-iextra_conv]

    #convolve
    date1ext = datesx[0] - timedelta(days=iextra_conv)
    datesx2 = list(rrule.rrule(rrule.DAILY, dtstart=date1ext, until=date2ext))

    x2 = np.array(datesx2)
    indx2= np.arange(x2.size)
    ftc = p * np.exp(k* (indx2-iextra_conv)) / (p * np.exp(k * (indx2-iextra_conv)) + (1 - p))
    minftc = 1-ftc

    ftc = lognormal_smooth_shift_convolve(ftc, gamma_mean, gamma_stdev, scale=1.0, fillval=ftc[0])
    ftc = lognormal_smooth_shift_convolve(ftc, -gamma_mean, 0, scale=1.0, fillval=ftc[0])

    minftc = lognormal_smooth_shift_convolve(minftc, gamma_mean, gamma_stdev, scale=1.0, fillval=ftc[0])
    minftc = lognormal_smooth_shift_convolve(minftc, -gamma_mean, 0, scale=1.0, fillval=ftc[0])

    ftc = ftc/(ftc+minftc)

    ftc = ftc[iextra_conv:-iextra_conv]

    return indx,inf,  ftc, Rt, r_lockdown, rtratio_voc, pu, pv, pb,  infpeak, hosp, hosppeak




def getRt(inf, ts):
    growth = np.diff(inf, prepend=inf[0])
    kgrowth = np.log( (growth+inf)/inf)
    gamma = 1/ts
    beta  = np.exp(kgrowth)-1+ gamma
    Rt = beta/gamma
    Rt[0]=Rt[1]
    return Rt

def getInterFt (ft, d, o):
    return (1-ft)*(d) + ft*(o)










def getpeaks(inf):
    infpeakdif = np.diff(inf)
    infstartgrowth = np.where(infpeakdif>0)
    if (np.size(infstartgrowth)>0):
        istart = infstartgrowth[0][0]
        infpeak = max (inf[istart:])
    else:
        infpeak = inf[-1]
    return infpeak





def base_sir_model(seir_init, Rt, ts):
    gamma = 1/ts
    s = seir_init[0]
    i = seir_init[1]
    r = seir_init[2]
    beta = Rt*gamma/seir_init[0]
    inf = np.empty([Rt.size])
    t = np.arange(inf.size)

    class BetaFunc(object):
        def __init__(self, beta, time):
            self.beta = beta
            self.time = time

        def get_beta(self, tx):
            # tx = min(tx,self.time[-1])
            # index = int (tx - self.time[0])
            # return self.alpha[index]
            return np.interp(tx, self.time, self.beta)

    beta_t = BetaFunc(beta, t)

    def ode(y, t,  gamma, ts):
        s,i, r = y

        # Basic SIR with time-dependent alpha

        dsdt = -beta_t.get_beta(t) * s * i
        didt =  beta_t.get_beta(t) * s * i - gamma * i
        drdt = gamma * i


        dydt = [dsdt,  # dS/dt      Susceptible
                didt,  # dI/dt      Infected
                drdt  # dR/dt      Removed
                ]
        #if (alpha > 0.98):
        #    print('alpha : ', alpha)
        return dydt

    init_vals = s,i,r
    tmin = t.min()
    tmax = t.max()
    lent = int(tmax - tmin)
    time_inflation = 1
    t_new = np.linspace(tmin, tmax, lent * time_inflation + 1)
    resultode = odeint(ode, init_vals, t_new,
                       args=(gamma, ts))
    infres = resultode[:,1]
    inf = np.interp(t,t_new, infres)
    return inf

def base_seir_model(seir_init, Rt,  k, ts, tlatent, tinf):
    gamma = 1.0/tinf
    sigma = 1.0/tlatent

    Rtref =Rt[0]

    # scale Rt such that it corresponds to the ts based Rt
    #lamba = Rt * (1/ts)- 1.0/ts
    #lamba = np.exp(k) - 1.0
    # Rseir = (lamba +sigma)*(lamba+gamma)/(sigma*gamma)
    #Rt *= Rseir/Rt
    #Rt -= Rt[0]-Rtref

    s = seir_init[0]
    i = seir_init[1]
    r = (1-s)
    beta = Rt * gamma/ s
    e = ((Rt[0]+1)*gamma/(2*sigma)) *i
    inf = np.empty([Rt.size])
    t = np.arange(inf.size)

    class BetaFunc(object):
        def __init__(self, beta, time):
            self.beta = beta
            self.time = time

        def get_beta(self, tx):
            # tx = min(tx,self.time[-1])
            # index = int (tx - self.time[0])
            # return self.alpha[index]
            return np.interp(tx, self.time, self.beta)

    beta_t = BetaFunc(beta, t)

    def ode(y, t, gamma, sigma):
        s, e, i, r = y
        dsdt =  -beta_t.get_beta(t)* s * i
        dedt =  beta_t.get_beta(t)* s * i - sigma * e
        didt = sigma * e - gamma * i
        drdt = gamma * i
        dydt = [dsdt,  # dS/dt      Susceptible
                dedt,  #dE/dt       Exposed
                didt,  # dI/dt      Infected
                drdt  # dR/dt      Removed
                ]
        return dydt

    init_vals = s, e, i, r
    tmin = t.min()
    tmax = t.max()
    lent = int(tmax - tmin)
    time_inflation = 1
    t_new = np.linspace(tmin, tmax, lent * time_inflation + 1)
    resultode = odeint(ode, init_vals, t_new,
                       args=(gamma, sigma))
    infres = resultode[:, 1]
    inf = np.interp(t, t_new, infres)

    return inf

def pdf_gamma( x ,shape,scale ):
    if (x<0):
        return 0
    else:
        return (x**(shape-1) *(np.exp(-x/scale)/ (sps.gamma(shape)*scale**shape)))

def pdf_lognormal (x, mu_log, sd_log):
    if (x>0):
        r = (1.0/(x*sd_log * np.sqrt(2*np.pi))) *  np.exp(- (np.log(x)- mu_log)**2 / (2*sd_log**2) )
    else:
        r= 0
    return r

def pdf_gamma( x ,shape,scale ):
    if (x<0):
        return 0
    else:
        return (x**(shape-1) *(np.exp(-x/scale)/ (sps.gamma(shape)*scale**shape)))



def find_lognormal_musd (mu_target, sd_target):
    mu = 0
    sd = 0
    mu_log = 0
    sd_log = 0.7


    while(abs(sd_target-sd)>1e-1):
        while (abs(mu_target -mu)>1e-2 ):
            mu = np.exp(mu_log + sd_log ** 2 / 2)
            mu_log = mu_log + np.log(mu_target/mu)
            sd = np.sqrt(np.exp(sd_log ** 2 + 2 * mu_log) * (np.exp(sd_log ** 2) - 1))
            #print ('mu, sd ', mu,sd)
        sd_log = sd_log + 0.3* np.log(sd_target /sd)
        mu = np.exp(mu_log + sd_log ** 2 / 2)
    return mu_log, sd_log



def lognormal_smooth_shift_convolve (input_nopad, mean, stddev, scale=1.0,fillval=0):
    """
    smooths the input with gaussian smooothing with standarddeviation and shifts its delay positions
    :param input: The input array
    :param mean: the amount of indices to shift the result
    :param the stddev for the gaussian smoothing (in index count)
    :param scale: scale the input array first with scale
    :return: the smoothed and shifted array
    """
    forcescale = False
    npad = 30
    input = np.pad(input_nopad, (npad,npad), 'edge')
    if isinstance(scale, np.ndarray):
        forcescale = True
    if (forcescale or np.abs(scale-1) > 1e-5):
        input = input*scale

    result = input
    if (stddev > 0.99) and (mean > 0.99):
        isd = 8* max(1,int (stddev))
        isd = min( int(0.5*np.size(input)-1), isd)
        mu_log, sd_log = find_lognormal_musd(mean, stddev)
        ishift = int (mean)
        mylognormal= np.fromiter(( pdf_lognormal(x, mu_log, sd_log) for x in range (-isd+ishift, isd+ishift+1)), np.float)
        result = np.convolve(input, mylognormal, mode='same')

    if (mean > 0):
        result = np.roll(result, int(mean))
        result[: int(mean)] = fillval
    else:
        result = np.roll(result, int(mean))
        # fill the trailing values with the last result
        result[mean:]=result[mean-1]

    result = result[npad:-npad]
    return result

def gamma_smooth_shift_convolve (input_nopad, mean, stddev, scale=1.0, fillval=0):
    """
    convolves the input with a gamma distribution with mean of shift
    :param input_nopad: The input array
    :param mean: the mean (choose integer to be sure that it is shifted the right way
    :param the stddev for the gaussian smoothing (in index count)
    :param scale: scale the input array first with scale
    :param shiftmean: shift the result with the mean
    :return: the convolved (and shifted) array
    """
    forcescale = False

    npad = 30
    input = np.pad(input_nopad, (npad,npad), 'edge')
    if isinstance(scale, np.ndarray):
        forcescale = True
    if (forcescale or np.abs(scale-1) > 1e-5):
        input = input*scale

    result = input
    if (stddev > 0.99) and (mean > 0.99):
        if (mean<stddev):
            stddev = mean - 1
        isd = 8* max(1,int (stddev))
        isd = min( int(0.5*np.size(input)-1), isd)
        theta = stddev ** 2 / mean
        k = mean / theta
        ishift = int (mean)
        mygamma= np.fromiter(( pdf_gamma(x, k, theta) for x in range (-isd+ishift, isd+ishift+1)), np.float)
        scale = np.sum (mygamma)
        mygamma /= scale
        result = np.convolve(input, mygamma, mode='same')

    if (mean > 0):
        result = np.roll(result, int(mean))
        result[: int(mean)] = fillval
    else:
        result = np.roll(result, int(mean))
        # fill the trailing values with the last result
        result[mean:]=result[mean-1]

    result = result[npad:-npad]
    return result




def  plot_ft(x,ft, perc, k):
    # Create figure and plot space
    fig, ax = plt.subplots(figsize=(20, 10))
    plt.xlabel('Date')
    plt.ylabel('VOC prevalence')
    ax.fill_between(x, ft[1], ft[2], facecolor='orange', label='range')
    slabel = "VOC prevalence f0 (median={},min={},max={}), k (median={},min={},max={}) ".format(perc[0], perc[1],
                                                                                                perc[2], k[0], k[1],
                                                                                                k[2])
    ax.plot(x,
            ft[0],
            color='green', label=slabel)
    plt.legend(loc='upper left')
    plt.grid()
    plt.show()


def plot_rt(x,Rt, Rtvoc, plotvac, agevac, demvac):
# Create figure of Rt evolution
    fig, ax = plt.subplots(figsize=(20, 10))
    plt.xlabel('Date')
    plt.ylabel('Rt')
    slabel = "Rratio (Rt(VOC)/Rt(ref))={}, Rt(ref)=Rlockdown={} ".format(Rtvoc[0],Rt_ld)
    ax.plot(x,
        Rt[0],
        color='mediumorchid', label=slabel)
    ax.fill_between(x, Rt[1], Rt[2], facecolor='plum')
    plt.legend(loc='center left')
    plt.grid()

    if (plotvac):
        ax2 = ax.twinx()
        ax2.plot(x,
            agevac,
            color='black', label='age vaccinated')
        percnotvac = 100- demvac *100/17.4
        ax2.plot(x,
            percnotvac,
            color='black', label='%not vaccinated', linestyle='dashed')
        plt.ylim(0,100)
        plt.ylabel('age or %')
        plt.legend(loc='lower right')
    plt.show()


def plot_inf(x, inf, hospsum, perc, plotvac,agevac,demvac):
    # create figure of growth of infection and relative hospitalized based on the agegroup differentation
    fig, ax = plt.subplots(figsize=(20, 10))
    plt.xlabel('Date')
    plt.ylabel('Ratio')

    for i, p in enumerate(perc):
        slabel = "infected Rlockdown={}, Rratio={}, f0={} ".format(Rt_ld, Rtvoc[0], p)
        if (i == 0):
            ax.plot(x,
                    inf[i],
                    color='lightcoral', label=slabel)
        else:
            ax.plot(x,
                    inf[i],
                    color='mistyrose', label=slabel, linestyle='dashed')
        slabel = "hospitalized Rlockdown={}, RtVoc={}, f0={} ".format(Rt_ld, Rtvoc[0], p)
        if (i == 0):
            ax.plot(x,
                    hospsum[i],
                    color='steelblue', label=slabel)
        else:
            ax.plot(x,
                    hospsum[i],
                    color='powderblue', label=slabel, linestyle='dashed')
    # ax.fill_between(x, hospsum[1], hospsum[2], facecolor='powderblue', label='conf', alpha=0.8)
    plt.legend(loc='center left')
    plt.ylim(0, 10)
    plt.grid()

    if (plotvac):
        ax2 = ax.twinx()
        ax2.plot(x,
                 agevac,
                 color='black', label='age vaccinated')
        percnotvac = 100 - demvac * 100 / 17.4
        ax2.plot(x,
                 percnotvac,
                 color='black', label='%not vaccinated', linestyle='dashed')
        plt.ylim(0, 100)
        plt.ylabel('age or %')
        plt.legend(loc='lower right')
    plt.show()



def single_run_method2_obsolete(r_lockdowndayx, r_lockdownval, r_lockdownscale,
               k_voc,  f0_voc, vac_ratio, startdate, enddate, boosterdayx, boosterdayn,
               p_vac_start,  p_booster_start,  demage, demn,
               agegroups_p_vac, w_hosp_age,
               startinfperc, startsusceptible,  ts,
               ve_trans = 0.5,
               ve_immune_o=0.2, ve_immune_hosp_o=0.5,
               ve_vac_d = 0.6, ve_vac_o=0.2,
               ve_booster_d=0.95, ve_booster_o=0.75,
               ve_vac_hosp_d= 0.95, ve_vac_hosp_o=0.48,
               ve_booster_hosp_d=0.975, ve_booster_hosp_o=0.9,
               re0_voc=-1):
    """
         construct Grid instance from fname and format

     :param  r_lockdowndayx: days after start specifying r_lockdown (for delta)
     :param  r_lockdownval: the R(0) value at r_lockdowndayx  (for delta)
     :param  r_lockdownscale  scaling value for the lockdown
     :param  k_voc : the daily growth rate (natural logarithm of the VOC (omicron))

     :param  vac_ratio: the  booster speed as fraction of boosterdayx, boosterdayn
     :param  boosterdayx: array with days from start
     :param  boosterdayn: number of boosters planned at boosterdayx
     :param  p_vac_start:  fraction of population vaccinated at start of simulation
     :param  p_booster_start:  fraction of population boostered at start of simulation

     :param  demage :  array with age cohort limits in the demography (descending, starting from ca 90 to 0)
     :param  demn : number of people in the age cohort
     :param  startinfperc: starting infection fraction of population (daily infected * ts)
     :param  startsusceptable: starting susceptable fraction of population (taking into account passed infections)
     :param  ts : serial interval time
     :param  Ud:  immune fraction in population related to
     :param  ve_trans  vaccin and booster efficy against transmission
     :param  ve_vac_d:  Vaccin efficacy against infections for delta
     :param  ve_vac_o   Vaccin efficacy against infections for omicron
     :param  ve_booster_d:  booster efficacy against infections for delta
     :param  ve_booster_o   booster efficacy against infections for omicron
     :param  re0_voc:  the ration of the basi reproduction number of the VOC omicron  to prevailing (delta), if omitted it will estimate it from the vaccination state and k value

     """

    date1 = startdate
    date2 = enddate
    datesx = list(rrule.rrule(rrule.DAILY, dtstart=parser.parse(date1), until=parser.parse(date2)))

    x = np.array(datesx)
    indx = np.arange(x.size)

    ft = np.empty([indx.size])


    k = k_voc







    p = 10**f0_voc

    ft = p * np.exp(k* indx) / (p * np.exp(k * indx) + (1 - p))


    demtotal = demn[-1]

    # number of boosters
    dembooster = np.interp(indx, boosterdayx, boosterdayn)
    dembooster *=vac_ratio

    # you cannot booster more than vaccinated
    dembooster = np.minimum(dembooster, demtotal*(p_vac_start))

    dp_b = dembooster/demtotal

    # unvaccinated kept fixed
    pu = ft*0 + (1.0 -  p_vac_start - p_booster_start)
    # create time series of pv and pb
    pv = ft * 0 + p_vac_start - dp_b
    pb = ft * 0 + p_booster_start + dp_b

    Vd = (1-ve_vac_d)*(1-ve_trans)
    Bd = (1-ve_booster_d)*(1-ve_trans)
    Vo = (1-ve_vac_o)*(1-ve_trans)
    Bo = (1-ve_booster_o)*(1-ve_trans)

    Ud = (1-startsusceptible)*(1-ve_immune_o) / startsusceptible


    #rtratio_voc, R0o_R0d = get_rtratio_voc(pu, pv, pb, Ud, Vo, Bo, Vd, Bd, p_vac_start, p_booster_start, k, ts)


    # eq 2 correction factors
    Cv =  (pu[0] + (p_vac_start + p_booster_start) * (1 - ve_trans))
    Ct_o = ft * 0 + ((pu + (pv + pb + Ud) * (1 - ve_trans)) / (pu + (pv + pb + Ud))) / Cv
    Ct_d = ft * 0 + (pu + (pv + pb) * (1 - ve_trans)) / Cv

    Id = ft * 0 + (pu + pv * (1 - ve_vac_d) + pb * (1 - ve_booster_d)) / (pu[0] + p_vac_start * (1 - ve_vac_d) + p_booster_start * (1 - ve_booster_d))
    Io = ft * 0 + (pu + (pv +Ud)* (1 - ve_vac_o) + pb * (1 - ve_booster_o)) / (pu[0] + (p_vac_start) * (1 - ve_vac_o) + p_booster_start * (1 - ve_booster_o))


    Rtd = Ct_d * Id

    gamma = 1 / ts
    beta = np.exp(k) - 1 + gamma
    rtratioref = beta / gamma
    R0o_R0d = rtratioref * Ct_d[0]*Id[0]/(Ct_o[0]*Io[0])
    #R0o_R0d = rtratioref * Ct_d[0]  / Ct_o[0]

    Rto=   R0o_R0d* Ct_o * Io
    rtratio_voc = Rto/Rtd

    rtratio_voc, R0o_R0d = get_rtratio_voc_obsolete(pu, pv, pb, Ud, Vo, Bo, Vd, Bd, p_vac_start, p_booster_start, k, ts)


    # apply effects of increase due to higher basic reproduction number
    Rtd = ft * 0 + 1
    Rtd = getInterFt(ft, Ct_d,  R0o_R0d *Ct_o)

    # now correct properly for R(0) and measures for additional lockdown
    r_lockdown = np.interp(indx, r_lockdowndayx, r_lockdownval)
    Rtd *= r_lockdown*r_lockdownscale/Rtd[0]

    # now apply future change in transmission due to shift towards boosters

    infected = startinfperc * 1e-2
    sirinit = [startsusceptible - infected, infected, 1- startsusceptible - infected]



    inf, pu,pv,pb = base_sirvb_model_obsolete(sirinit, ve_immune_o, Rtd,   ts, ft, pu, pv, pb, ve_vac_d, ve_vac_o, ve_booster_d,   ve_booster_o,
                                     ve_trans)



    hosp = getHosp(inf, ft, demage, demn, agegroups_p_vac, w_hosp_age, p_booster_start, p_vac_start, dp_b,
            Ud, ve_immune_hosp_o, ve_vac_hosp_d, ve_vac_hosp_o, ve_booster_hosp_d, ve_booster_hosp_o)

    inf = inf / inf[0]
    infpeak = getpeaks(inf)

    Rt = getRt(inf, ts)


    # obtain the the hospitalized fraction as a function of different vacinne status and




    hosppeak = getpeaks(hosp)

    return indx,inf,  ft, Rt, Rtd, rtratio_voc, pu, pv, pb,  infpeak, hosp, hosppeak


def get_rtratio_voc_obsolete(pu, pv, pb, Ud, Vo, Bo, Vd, Bd, p_vac_start, p_booster_start, k, ts):
    rtratio_voc = (pu + (Ud + pv) * Vo + pb * Bo) / (pu + pv * Vd + pb * Bd)

    gamma = 1 / ts
    beta = np.exp(k) - 1 + gamma
    rtratioref = beta / gamma

    rtratiovacstate = (pu + (Ud + p_vac_start) * Vo + p_booster_start * Bo) / (
                pu[0] + p_vac_start * Vd + p_booster_start * Bd)

    R0o_R0d = rtratioref / rtratiovacstate

    rtratio_voc *= R0o_R0d

    #print (" R0o_R0d ",R0o_R0d[0] )

    return rtratio_voc, R0o_R0d



def base_sirvb_model_obsolete(seir_init, ve_immune_o, Rtref,  ts, ft, pu, pv, pb, vinf_d, vinf_o, binf_d, binf_o, ve_trans):
    gamma = 1/ts
    s = s0 = seir_init[0]*pu[0]
    i = seir_init[1]
    r = seir_init[2]
    v = v0=  seir_init[0]*pv[0]
    b = b0 = seir_init[0]*pb[0]
    eta_v = getInterFt(ft[0], 1.0-vinf_d, 1.0-vinf_o)
    eta_b = getInterFt(ft[0], 1.0-binf_d, 1.0-binf_o)
    svb = s+ v*eta_v + b*eta_b
    beta = Rtref*gamma/svb

    inf = np.empty([Rtref.size])
    t = np.arange(inf.size)

    class BetaFunc(object):
        def __init__(self, beta, time):
            self.beta = beta
            self.time = time

        def get_beta(self, tx):
            return np.interp(tx, self.time, self.beta)

    eta_vv = getInterFt(ft, 1.0 - vinf_d, 1.0 - vinf_o)
    eta_bb = getInterFt(ft, 1.0 - binf_d, 1.0 - binf_o)



    beta_t = BetaFunc(beta, t)
    eta_v_t = BetaFunc(eta_vv, t)
    eta_b_t = BetaFunc(eta_bb, t)
    dfdt_t = BetaFunc(np.diff(ft), t[0:-1])
    dvdt_t = BetaFunc(np.diff(pv), t[0:-1])
    dbdt_t = BetaFunc(np.diff(pb), t[0:-1])

    def ode(y, t,  gamma, ts):
        s, i, r, v, b = y

        # Basic SIR with time-dependent alpha

        bb = beta_t.get_beta(t)  # * ((s + (v +b)*(1-ve_trans))*(s0+v0+b0))/((s0+ (v0+b0)*(1-ve_trans))*(s+b+v))
        dsdt = -bb * s * i
        didt = bb *  i * (s + v * eta_v_t.get_beta(t)  + b  * eta_b_t.get_beta(t)) - gamma * i
        dvdt = dvdt_t.get_beta(t)*(1-r) + dfdt_t.get_beta(t)* (seir_init[2]*(1-ve_immune_o)) - bb *  i * v * eta_v_t.get_beta(t)
        dbdt = dbdt_t.get_beta(t)*(1-r)  - bb *  i * b * eta_b_t.get_beta(t)
        drdt = gamma * i - dfdt_t.get_beta(t)* (seir_init[2]*(1-ve_immune_o))


        dydt = [dsdt,  # dS/dt      Susceptible
                didt,  # dI/dt      Infected
                drdt,  # dR/dt      Removed
                dvdt,
                dbdt
                ]
        #if (alpha > 0.98):
        #    print('alpha : ', alpha)
        return dydt

    init_vals = s,i,r,v,b
    tmin = t.min()
    tmax = t.max()
    lent = int(tmax - tmin)
    time_inflation = 1
    t_new = np.linspace(tmin, tmax, lent * time_inflation + 1)
    resultode = odeint(ode, init_vals, t_new,
                       args=(gamma, ts))
    infres = resultode[:,1]
    rr = resultode[:, 2]
    ss = resultode[:, 0]
    vv = resultode[:, 3]
    bb = resultode[:, 4]
    inf = np.interp(t,t_new, infres)
    sum = np.interp(t,t_new, ss) + np.interp(t,t_new, vv)+ np.interp(t,t_new, bb)
    pu = np.interp(t,t_new, ss) / sum
    pv = np.interp(t,t_new, vv) / sum
    pb = np.interp(t,t_new, bb) / sum
    return inf,  pu, pv, pb



