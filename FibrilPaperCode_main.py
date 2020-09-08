'''
Created on Jan 6, 2020

@author: juliaziaee
'''

import csv
from numpy import double
import numpy as np
from scipy import stats
from sklearn.metrics import r2_score

'''
reads in csv file to get rates
'''


def filereadrates(filename):
    rates = []
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        fields = next(csvreader)
        for row in csvreader:
            if row[0] != '':
                rates.append(row)

    return rates


'''
reads in csv file to get time points
'''


def filereadtimes(filename):
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        fields = next(csvreader)
        times = (fields)
    times[0] = '0'

    return times


'''
sets up output using the other functions in the program
'''


def calc(filename):
    times = filereadtimes(filename)
    rates = filereadrates(filename)

    cutoffs = ratechecker(rates)
    opt = optimalRate(rates)
    ct = ctvals(rates, cutoffs, times)
    maxtime = times[opt]
    ct_by_med = manipulation_median(ct[1])
    ct_by_mean = manipulation_mean(ct[1])
    return "Threshold Value is " + str(ct[0]) + " and optimal linearity is at time " + str(
        maxtime) + " \n ct vals in order: " + str(ct[1]) + "\n median ct value for each sample: " + str(
        ct_by_med) + "\n mean ct value for each sample: " + str(ct_by_mean) + "\n R-Square for the means: " + str(ct[2])


'''
Finds threshold value used to get Ct values
'''


def ratechecker(rates):
    triplets = []
    strtind = 0
    builder = []
    for rate in rates:
        if strtind % 3 != 0 or strtind == 0:
            builder.append(rate)
            strtind += 1
        else:
            triplets.append(builder)
            builder = []
            builder.append(rate)
            strtind += 1

    hitnum = 0.1
    thresholds_to_check = []
    while hitnum <= 0.3:
        group_ct = []
        for trip in triplets:
            avgs = []
            for t in trip:
                i = 1
                comps = []
                while i < len(t):
                    top = double(t[i])
                    bottom = double(t[i - 1])
                    rtchange = top / bottom
                    comps.append(rtchange)
                    i += 1
                n = 1
                checker = 0
                while n < len(comps) - 1 and checker == 0:
                    percentchange = comps[n] - comps[n - 1]
                    if percentchange >= hitnum and comps[n] > 1 and comps[n + 1] > 1:
                        avgs.append(double(t[n + 1]))
                        checker += 1
                    n += 1

            look = sorted(avgs)
            if len(look) == 3:
                group_ct.append(look[1])
            if len(look) == 2:
                group_ct.append((look[0] + look[1]) / 2)
            if len(look) == 1:
                group_ct.append(look[0])
            if len(look) == 0:
                samp = trip[0]
                sampnum = samp[len(samp) - 1]
                group_ct.append(double(sampnum))

        helpersum = 0
        if len(group_ct) != 0:
            for num in group_ct:
                helpersum += num

        hitnum += 0.01
        thresholds_to_check.append(helpersum / len(group_ct))

    return thresholds_to_check


'''
Finds time of optimal linearity
'''


def optimalRate(rates):
    wininds = []
    for rate in rates:
        rtchanges = []
        correspond_ind = []
        strtind = 1
        while strtind < len(rate):
            change = double(rate[strtind]) - double(rate[strtind - 1])
            rtchanges.append(change)
            correspond_ind.append(strtind)
            strtind += 1
        winnerind = 0
        winnerrt = 0
        for i in range(len(rtchanges)):
            if rtchanges[i] > winnerrt:
                winnerrt = rtchanges[i]
                winnerind = correspond_ind[i]
        wininds.append(winnerind)
    result = sorted(wininds, key=wininds.count, reverse=True)
    return result[0]


'''
Finds ct values for each triplicate in the csv
'''


def ctvals(rates, threshs, times):
    winning_thresh = 0
    winning_ct_vals = []
    winning_rsq = 0

    for thresh in threshs:
        siind = 0
        t = []
        intermediate = []
        for rate in rates:
            if len(intermediate) == 3:
                t.append(intermediate)
                intermediate = []
            ind = 0
            done = 0
            for num in rate:
                if double(num) >= thresh and done == 0 and ind != 0:
                    intermediate.append(double(times[ind]))
                    done += 1
                    ind += 1
                else:
                    if ind == (len(rate) - 1) and done == 0:
                        intermediate.append(double(times[(len(rate) - 1)]))
                    ind += 1
            if siind == len(rates) - 1:
                t.append(intermediate)
            siind += 1

        means_to_eval = manipulation_mean(t)
        rval = rsquare(means_to_eval)
        if rval > winning_rsq:
            winning_rsq = rval
            winning_ct_vals = t
            winning_thresh = thresh

    return [winning_thresh, winning_ct_vals, winning_rsq]


'''
returns median of the ct vals found for each sample
'''


def manipulation_median(trips):
    med = []
    for trip in trips:
        trip.sort()
        med.append(trip[1])
    return med


'''
returns mean of the ct vals found for each sample
'''


def manipulation_mean(trips):
    men = []
    for trip in trips:
        s = 0
        for num in trip:
            s += num
        avgval = round(s / len(trip), 2)
        men.append(avgval)
    return men


def rsquare(means):
    xs = []
    for i in range(len(means)):
        xs.append(i)
    # mean = sorted(means)
    y = np.array(means)
    x = np.array(xs)

    r2 = r2_score(y, linefitline(x, x, y))

    return r2


def linefitline(b, x, y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return (intercept + slope * b)


filename = "Arpi_SD.csv"
print(calc(filename))