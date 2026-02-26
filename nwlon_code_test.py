# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 11:00:00 2026

@author: David Ilogho
"""

from sl3 import *
import serial
from _thread import allocate_lock as locker
from uos import listdir
from urandom import uniform
import math

# PARAMETERS — customize as needed
TIDE_PERIOD_SECONDS = 12 * 60 * 60
TIDE_AMPLITUDE = 3
SINE_NOISE_FREQ = 0.25
SINE_NOISE_AMP  = 1.5


def live_data(add_constant=15):
    curr_time = int(utime.time())
    cycle_pos = (curr_time % TIDE_PERIOD_SECONDS) / TIDE_PERIOD_SECONDS
    angle = 2 * math.pi * cycle_pos
    base_tide = TIDE_AMPLITUDE * math.sin(angle)
    noise = SINE_NOISE_AMP * math.sin(2 * math.pi * SINE_NOISE_FREQ * curr_time)
    if uniform(0,1) > 0.97:
        return  base_tide + noise + add_constant
    return base_tide + noise + uniform(0,5)

# Declaring variables
add_sns = []
ag_active = aq_active = dqap1_init = dqap2_init = False
cnt_meas = curr_time = 0
dqap1_class = dqap2_class = m1_sns_check = paros_class = vis_class = None
goes_msg = ''
gp5_status_msg = float(command_line('!GP5 Value\r'))
initialize = paros_init = vis_init = False
lock_write = locker()
Log.limit = 1000000
no_of_meas = 32
station_id = command_line('!STATION NAME\r').strip()
six_min_avg = 'MWWL', 'BWL', 'AG'
three_min_avg = 'AQT', 'PAROS', 'SAE', 'VIS3M'
two_min_avg = 'AT', 'BARO', 'COND', 'RH', 'WT'

aq_temp = (
    ['AQT', 'AQTSTD', 'AQTOUT'],
    ['AQT1', 'AQT2']
)

ag_temp = (
    ['AG', 'AGSTD', 'AGOUT'],
    ['AGT1', 'AGT2']
)

aq = aq_temp[0]
ag = ag_temp[0]

all_sns = (
    aq,
    ('SAE', 'SAESTD', 'SAEOUT'),
    ('SAE2', 'SAESTD2', 'SAEOUT2'),
    ('MWWL', 'MWSTD', 'MWOUT'),
    ('MWWL2', 'MWSTD2', 'MWOUT2'),
    ('PAROS', 'PAROSSTD', 'PAROSOUT'),
    ('PAROS2', 'PAROSSTD2', 'PAROSOUT2'),
    ('BPAROS', 'BPAROSSTD', 'BPAROSOUT'),
    ('BPAROS2', 'BPAROSSTD2', 'BPAROSOUT2'),
    ag,
    ('AG2', 'AGSTD2', 'AGOUT2'),
    ('WS', 'WD', 'WG'),
    ('WS2', 'WD2', 'WG2'),
    ('AT',), ('AT2',),
    ('WT',), ('WT2',),
    ('BARO',), ('BARO2',),
    ('COND',), ('COND2',),
    ('RH',), ('RH2',),
    ('VIS3M', 'ALV'),
    ('AUXBAT',), ('BAT',),
    ('BWL', 'BWLSTD', 'BWLOUT'),
    ('BBAT',),
    ('SNS', 'DAT'),
    ('AQTWL', 'MWTWL')
)

all_sns_remaining = (('SNS',), ('DAT',), ('AQTWL',), ('MWTWL',))

all_sns_single = (
    'AQT', 'AQTSTD', 'AQTOUT', 'AQT1', 'AQT2',
    'SAE', 'SAESTD', 'SAEOUT',
    'SAE2', 'SAESTD2', 'SAEOUT2',
    'MWWL', 'MWSTD', 'MWOUT',
    'MWWL2', 'MWSTD2', 'MWOUT2',
    'PAROS', 'PAROSSTD', 'PAROSOUT',
    'PAROS2', 'PAROSSTD2', 'PAROSOUT2',
    'BPAROS', 'BPAROSSTD', 'BPAROSOUT',
    'BPAROS2', 'BPAROSSTD2', 'BPAROSOUT2',
    'AG', 'AGSTD', 'AGOUT', 'AGT1', 'AGT2',
    'AG2', 'AGSTD2', 'AGOUT2',
    'WS', 'WD', 'WG', 'WS2', 'WD2', 'WG2',
    'AT', 'AT2', 'WT', 'WT2',
    'BARO', 'BARO2', 'COND', 'COND2',
    'RH', 'RH2', 'VIS3M', 'ALV',
    'AUXBAT', 'BAT', 'BBAT',
    'BWL', 'BWLSTD', 'BWLOUT',
    'SNS', 'DAT',
    'AQTWL', 'MWTWL'
)

baro = 'BARO', 'BARO2'

discard_sns = (
    'AQTCOUNTS',
    'BWLCOUNTS',
    'DQAP',
    'SAECOUNTS', 'SAECOUNTS2',
    'MWCOUNTS', 'MWCOUNTS2',
    'PCOUNTS', 'PCOUNTS2',
    'BPCOUNTS', 'BPCOUNTS2',
    'AGCOUNTS', 'AGCOUNTS2',
    'PAROSRAW',
    'VIS10M', 'VISRAW',
    'WndSpd', 'WndDir',
    'WndSpd2', 'WndDir2'
)

m1_sns = (
    'AQT', 'SAE', 'MWWL',
    'PAROS', 'BPAROS', 'AG',
    'WS', 'AT', 'WT', 'BARO',
    'COND', 'RH', 'VIS3M', 'BWL'
)

must_include_sns = 'SNS', 'DAT', 'BAT'

one_byte_sns = (
    'AQTOUT',
    'SAEOUT', 'SAEOUT2',
    'MWOUT', 'MWOUT2',
    'BWLOUT',
    'PAROSOUT', 'PAROSOUT2',
    'BPAROSOUT', 'BPAROSOUT2',
    'AGOUT', 'AGOUT2',
    'ALV'
)

primary_sns = (
    'AQT',
    'BWL',
    'SAE', 'SAE2',
    'MWWL', 'MWWL2',
    'PAROS', 'PAROS2',
    'BPAROS', 'BPAROS2',
    'AG', 'AG2'
)

primary_sns_grp = (
    aq,
    ('BWL', 'BWLSTD', 'BWLOUT'),
    ('WS', 'WD', 'WG'),
    ('SAE', 'SAESTD', 'SAEOUT'),
    ('MWWL', 'MWSTD', 'MWOUT'),
    ('PAROS', 'PAROSSTD', 'PAROSOUT'),
    ('BPAROS', 'BPAROSSTD', 'BPAROSOUT'),
    ag,
    ('VIS3M', 'ALV')
)

primary_sns_grp_remaining = (
    aq,
    ('BWL', 'BWLSTD', 'BWLOUT', 'BBAT'),
    ('SAE2', 'SAESTD2', 'SAEOUT2'),
    ('MWWL2', 'MWSTD2', 'MWOUT2'),
    ('PAROS2', 'PAROSSTD2', 'PAROSOUT2'),
    ('BPAROS2', 'BPAROSSTD2', 'BPAROSOUT2'),
    ('AG2', 'AGSTD2', 'AGOUT2'),
    ('WS2', 'WD2', 'WG2')
)

tsunami_sns = 'AQTWL', 'MWTWL'

two_byte_sns = (
    'AQTSTD',
    'SAESTD', 'SAESTD2',
    'MWSTD', 'MWSTD2',
    'PAROSSTD', 'PAROSSTD2',
    'BPAROSSTD', 'BPAROSSTD2',
    'AGSTD', 'AGSTD2',
    'WS', 'WD', 'WG', 'WS2', 'WD2', 'WG2',
    'BARO', 'BARO2',
    'BWLSTD',
    'BAT', 'AUXBAT', 'BBAT',
    'RH', 'RH2'
)

two_byte_sns_sign = (
    'AT', 'AT2', 'WT', 'WT2', 'AQT1', 'AQT2', 'AGT1', 'AGT2', 'SNS'
)


class SecondarySensor:
    """
    This SecondarySensor class constructor uses the measurement number e.g.
    M1 to retrieve relevant sensor attributes, then uses it to create and
    initialize the sensor object.
    :param meas_number: Measurement number e.g. M1, M2.
    """

    def __init__(self, meas_number):
        """
        This SecondarySensor class constructor initializes the
        secondary sensor parameters with default values.
        """
        self.meas_number = meas_number
        self.label = setup_read('{0} Label'.format(meas_number))
        self.right_digits = int(
            setup_read('{0} Right digits'.format(meas_number))
        )
        self.value = -99999.0

    def get_encoded_data(self):
        """
        This method returns the encoded data of the object when called.
        :return: Encoded data
        """
        if self.value == -99999.0:
            value = 999999.0
        else:
            value = int(self.value * 10 ** self.right_digits)
        if self.label in one_byte_sns:
            return pseudo_encoder(value, 1)
        if self.label in two_byte_sns:
            if self.label in baro:
                value -= 8000
            return pseudo_encoder(value, 2)
        if self.label in two_byte_sns_sign:
            return pseudo_encoder(value, 2, True)
        return pseudo_encoder(value, 3)

    def update_secondary_data(self):
        """
        This method is used to update the secondary data object
        with the most recent sensor data.
        :return: Most recent sensor data value.
        """
        old = curr_time - 360
        new = curr_time
        log = [r for r in Log(count=1, match=self.label, oldest=old, newest=new)]
        if len(log) == 1:
            self.value = round(log[0].value, self.right_digits)
        else:
            self.value = -99999.0
        return self.value


class PrimarySensor(SecondarySensor):
    """
    This PrimarySensor class constructor inherits the SecondarySensor class
    attributes and methods. It also includes date and time attributes and
    additional methods for date and time conversions.
    :param meas_number: Measurement number e.g. M1.
    """

    def __init__(self, meas_number):
        """
        This PrimarySensor class constructor initializes the
        primary sensor parameters with default values.
        """
        super().__init__(meas_number)
        self.redundant_value = -99999.0
        init_time = utime.localtime()
        self.year, self.month, self.day = init_time[0], init_time[1], init_time[2]
        self.hour, self.minute, self.second = init_time[3], init_time[4], init_time[5]
        self.julian_day = init_time[7]
        self.sutron_day = sutron_day_calc(self.julian_day, self.year)

    def get_encoded_hour(self):
        """
        This method encodes hour.
        :return: returns encoded data.
        """
        return pseudo_encoder(self.hour, 1)

    def get_encoded_minute(self):
        """
        This method encodes minute.
        :return: returns encoded data.
        """
        return pseudo_encoder(self.minute, 1)

    def get_encoded_redundant_data(self):
        """
        This method encodes redundant data.
        :return: returns encoded data.
        """
        if self.redundant_value == -99999.0:
            value = 999999.0
        else:
            value = int(self.redundant_value * 10 ** self.right_digits)
        return pseudo_encoder(value, 3)

    def get_encoded_sutron_day(self):
        """
        This method encodes sutron day.
        :return: returns encoded data.
        """
        return pseudo_encoder(self.sutron_day, 2)

    def update_primary_data(self):
        """
        This method is used to update the primary data object with the most
        recent sensor data, date and time.
        :return: Most recent sensor data value, date and time.
        """

        # Retrieve most recent primary data
        old = curr_time - 360
        new = curr_time
        log = [r for r in Log(count=1, match=self.label, oldest=old, newest=new)]
        if len(log) == 1:
            self.value = round(log[0].value, self.right_digits)

            if self.label in six_min_avg:
                time_diff = -180
            elif self.label in three_min_avg:
                time_diff = -90
            elif self.label in two_min_avg:
                time_diff = -60
            else:
                time_diff = 0
            pri_date = utime.localtime(log[0].time + time_diff)
            self.year = pri_date[0]
            self.month = pri_date[1]
            self.day = pri_date[2]
            self.hour = pri_date[3]
            self.minute = pri_date[4]
            self.second = pri_date[5]
            self.julian_day = pri_date[7]
            self.sutron_day = sutron_day_calc(self.julian_day, self.year)

            # Retrieve redundant primary data
            old -= 360
            new -= 360
            log = [r for r in Log(count=1, match=self.label, oldest=old, newest=new)]
            if len(log) == 1:
                self.redundant_value = round(log[0].value, self.right_digits)
            else:
                self.redundant_value = -99999.0
        else:
            self.value = self.redundant_value = -99999.0


class TsunamiData(SecondarySensor):
    """
    This TsunamiData class constructor inherits certain PrimarySensor class
    attributes and methods. It also includes attributes and methods that are
    tsunami related.
    :param meas_number: Measurement number e.g. M1.
    """

    def __init__(self, meas_number):
        """
        This Tsunami class constructor initializes the tsunami sensor parameters
        with default values.
        """
        super().__init__(meas_number)
        self.value2 = self.value3 = self.value4 = self.value5 = self.value6 = self.value
        init_time = utime.localtime()
        self.hour, self.minute = init_time[3], init_time[4]

    def get_encoded_tsunami(self):
        """
        This method encodes tsunami data.
        :return: returns encoded data.
        """

        def tsu_calc(val, rnd=False):
            """
            This function helps with the Tsunami encoding algorithm calculation.
            rnd stands for rounding.
            """
            if rnd:
                return round(val * 10 ** self.right_digits)
            return pseudo_encoder(
                val % 250 + (val // 250 - tsu_offset) * 250, 2
            )
        tsu_val = (
                self.value,
                self.value2,
                self.value3,
                self.value4,
                self.value5,
                self.value6
        )

        if -99999.0 in tsu_val:
            return '?', '?', '@', '@@', '@@', '@@', '@@', '@@', '@@'
        tsu_val = [tsu_calc(tv, True) for tv in tsu_val]
        tsu_offset = min(tsu_val) // 250
        tsu_val = [tsu_calc(tv) for tv in tsu_val]
        tsu_hour = pseudo_encoder(self.hour, 1)
        tsu_minute = pseudo_encoder(self.minute, 1)
        tsu_offset = pseudo_encoder(tsu_offset, 1)
        return (tsu_hour, tsu_minute, tsu_offset) + tuple(tsu_val)

    def update_tsunami_data(self):
        """
        This method is used to update the tsunami data object with the most
        recent tsunami sensor data, date and time.
        :return: Most recent tsunami sensor data value, date and time.
        """
        old = curr_time - 60
        new = curr_time
        val_temp = [-99999.0] * 6
        log = []
        for vt in range(6):
            log.append([r for r in Log(count=1, match=self.label, oldest=old, newest=new)])
            old -= 60
            new -= 60
            if len(log[-1]) == 1:
                val_temp[vt] = round(log[-1][0].value, self.right_digits)

        if len(log[0]) == 1:
            tsu_date = utime.localtime(log[0][0].time)
            self.hour = tsu_date[3]
            self.minute = tsu_date[4]

        (
            self.value,
            self.value2,
            self.value3,
            self.value4,
            self.value5,
            self.value6
        ) = val_temp


class DQAP:
    """
    This class calculates DQAP for logged measurements.
    """

    def __init__(self, log_duration, obj_name):
        """
        This DQAP class constructor initializes the DQAP parameters
        with default values.
        """
        self.log_duration = min(1800, log_duration)
        self.mean = self.std = self.out = self.count = -99999.0
        self.dqap_update = self.mean_update = 0.0
        self.object_name = obj_name
        self.old_meas_time = utime.time()
        self.dqap_list = []
        self.dqap_lock = locker()
        self.vis3 = self.vis10 = self.alv = -99999.0

    def build_list(self, data_recent, meas_obj):
        """
        This method is used to build the list for DQAP and visibility
        processing.
        """
        self.dqap_lock.acquire()
        try:
            self.dqap_update = meas_as_reading(meas_obj).time
            self.dqap_list.append((self.dqap_update, data_recent))
            if self.object_name == 'VIS':
                self.dqap_list[:] = self.dqap_list[-5:]
            else:
                self.dqap_list[:] = self.dqap_list[-(self.log_duration + 25):]
        finally:
            self.dqap_lock.release()

    def dqap_parameter(self, para_num):
        """
        This method is what updates visibility parameters.
        """
        deadline = utime.time() + 15
        para_name = 'standard deviation', 'outliers', 'counts'
        while utime.time() < deadline:
            if utime.time() - self.mean_update < 5:
                para_val = self.std, self.out, self.count
                if para_val[para_num] == -99999.0:
                    status_message(
                        '{0} {1} update unsuccessful!'.format(
                            self.object_name, para_name[para_num]
                        )
                    )
                else:
                    status_message(
                        '{0} {1} update successful!'.format(
                            self.object_name, para_name[para_num]
                        )
                    )
                return para_val[para_num]
            utime.sleep_ms(100)
        status_message(
            '{0} {1} update unsuccessful, waited too long...'.format(
                self.object_name, para_name[para_num]
            )
        )
        return -99999.0

    def reset(self):
        """
        This method initializes DQAP class.
        """
        if self.object_name == 'VIS':
            self.vis3 = self.vis10 = self.alv = -99999.0
        else:
            self.mean = self.std = self.out = self.count = -99999.0
        self.mean_update = self.old_meas_time = utime.time()

    def tsunami_avg(self, meas_obj):
        """
        This function checks for updated raw values and returns the updated
        values.
        """
        return update_dqap(self, 'TWL', meas_obj, 12, 59)

    def update_dqap(self, meas_obj):
        """
        This method checks for updated raw values and returns the updated
        values.
        """
        update_dqap(self, self.object_name, meas_obj, 12, self.log_duration)
        return self.mean


class VisData(DQAP):
    """
    This class inherits the DQAP class and utilizes some of its attributes
    and methods to process vis data.
    """

    def __init__(self):
        super().__init__(log_duration=150, obj_name='VIS')

    def return_vis_status(self):
        """
        This method returns saves the visibility status in a file.
        """
        deadline = utime.time() + 15
        while utime.time() < deadline:
            if self.dqap_lock.acquire(False):
                try:
                    status_message(return_vis_status())
                    return
                finally:
                    self.dqap_lock.release()
            utime.sleep_ms(50)
        status_message(
            'Visibility status update unsuccessful. Please try again!'
        )

    def update_vis(self, meas_obj):
        """
        This method retrieves the data from the list and updates visibility
        parameters.
        """
        update_dqap(self, self.object_name, meas_obj, 10, self.log_duration)
        return self.vis3

    def vis_parameter(self, para_num):
        """
        This method retrieves already processed visibility data.
        """
        deadline = utime.time() + 13
        para_name = 'VIS10M', 'ALV'
        while utime.time() < deadline:
            if utime.time() - self.mean_update < 5:
                para_val = self.vis10, self.alv
                if para_val[para_num] == -99999.0:
                    status_message(
                        '{0} update unsuccessful...'.format(
                            para_name[para_num]
                        )
                    )
                else:
                    status_message(
                        '{0} update successful!'.format(para_name[para_num])
                    )
                return para_val[para_num]
            utime.sleep_ms(100)
        status_message(
            '{0} update not successful!'.format(para_name[para_num])
        )
        return -99999.0


def calc_dqap(data, only_mean=False, k=3.0):
    """
    Return [mean, stdev, outliers, counts] using k-sigma rule.
    """
    vals = [float(x) for x in data]
    total = len(vals)
    if total < 2:
        if only_mean:
            return -99999.0
        return [-99999.0] * 4

    # raw mean & stdev
    mu = sum(vals) / total
    if only_mean:
        return mu
    var = sum((x - mu) ** 2 for x in vals) / (total - 1)
    sd = var ** 0.5

    # remove outliers
    lo, hi = mu - k * sd, mu + k * sd
    cleaned = [x for x in vals if lo <= x <= hi]
    if not cleaned:  # if everything filtered out, keep original
        cleaned = vals[:]

    # recompute mean & stdev on cleaned set
    n = len(cleaned)
    dqap_mean = sum(cleaned) / n
    var = sum(
        (x - dqap_mean) ** 2 for x in cleaned
    ) / (n - 1) if n > 1 else 0.0
    sd = var ** 0.5
    return [dqap_mean, sd, total - n, n]


def calc_mean(data_list):
    """
    This function calculates and returns mean.
    """
    return calc_dqap(data_list, True)


def calc_vis(data_list):
    """
    This function returns visibility data.
    """
    data_list[:] = data_list[-3:]
    if len(data_list) == 0:
        return (-99999.0,) * 3
    viserr, vis3m, vis10m = zip(*data_list)
    vis3m = [v for v in vis3m if v >= 0]
    n = len(vis3m)
    if n > 0:
        return sum(vis3m) / n, vis10m[-1], viserr[-1]
    return -99999.0, vis10m[-1], viserr[-1]


def clear_dir():
    """
    This function clears the satlink 3 directory.
    """
    files_to_del = listdir('/')
    for ftd in files_to_del:
        try:
            if ftd.upper() not in ('SCRIPT.PY', 'SD'):
                command_line('!file del {0}\r'.format(ftd))
        except Exception:
            pass
    command_line('!script save\r')


def check_label_err(temp_label):
    """
    This function checks for incorrectly labeled measurements.
    Will not initialize if all condition are not met.
    """
    missing = []
    typo_label = []
    must_include = []
    err_config = False

    def p3_feedback():
        """
        This function appends configuration errors
        to p3 and updates error status.
        """
        nonlocal err_config
        err_config = True
        status_message(sm)
        with open('p3', 'a') as f:
            f.write('{0}: {1}\r\n'.format(sl3_datetime(), sm))

    # Check if measurement M1 is enabled
    if m1_sns_check is None:
        sm = 'Please ensure measurement M1 is enabled and reboot...'
        p3_feedback()

    # Check if measurement M1 has the correct measurement label
    elif m1_sns_check not in m1_sns:
        sm = 'Please ensure M1 has one of the following labels and reboot:\r\n  {0}'.format(
            m1_sns
        )
        p3_feedback()

    # Check if there are labels that aren't part of the accepted labels
    for tl in temp_label:
        if tl not in all_sns_single:
            typo_label.append(tl)

    if len(typo_label) > 0:
        for tl in typo_label:
            sm = 'Please verify {0} label and reboot, will not initialize at this time...'.format(
                tl
            )
            p3_feedback()

    # Check for incomplete group labels
    for psg in primary_sns_grp[2:] + primary_sns_grp_remaining:
        gl_length = len(psg)
        count = 0
        for p in psg:
            if p in temp_label:
                count += 1
        if count == 0 or gl_length == count:
            pass
        else:
            missing.append(psg)

    if len(missing) > 0:
        for m in missing:
            sm = 'Please verify {0} labels and reboot, will not initialize at this time...'.format(
                tuple(m)
            )
            p3_feedback()

    # Check for labels that must be included
    for mis in must_include_sns:
        if mis not in temp_label:
            must_include.append(mis)

    if len(must_include) > 0:
        for mi in must_include:
            sm = 'Please add {0} label and reboot, will not initialize at this time...'.format(
                mi
            )
            p3_feedback()

    return err_config


def init_paros():
    """
    This function executes the paros configuration
    if the paros sensor is enabled.
    """
    status_message('Configuring Paros sensor...')
    paros_reply = paros_baudset_9600()
    paros_reply = paros_setup(paros_reply)
    for pr in paros_reply[1]:
        status_message(pr)
    if paros_reply[0]:
        status_message(
            'Paros sensor configuration successful!'
        )
    else:
        status_message(
            'Paros sensor configuration was not successful!'
        )


def init_vis():
    """
    This function executes the visibility configuration
    if the visibility sensor is enabled.
    """
    status_message('Configuring Visibility sensor...')
    vis_reply = vis_initialize()
    for vr in vis_reply[1]:
        status_message(vr)

    if vis_reply[0]:
        status_message(
            'Visibility sensor configuration successful!'
        )
        vis_msg = '{0}: {1}\r\n'.format(
            sl3_datetime(), 'Visibility sensor configuration successful!'
        )
        with open('vs', 'w') as f:
            f.write(vis_msg)
    else:
        status_message(
            'Visibility sensor configuration not successful!'
        )
        vis_msg = '{0}: {1}\r\n'.format(
            sl3_datetime(), 'Visibility sensor configuration not successful...'
        )
        with open('vs', 'w') as f:
            f.write(vis_msg)


def paros_baudset_9600():
    """
    This function is used for setting paros baudrate.
    """
    bauds = 115200, 57600, 38400, 19200, 9600, 4800, 2400, 1200, 600, 300
    tm_out = 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2
    baud = tm = None
    for baud, tm in zip(bauds, tm_out):
        try:
            with serial.Serial('RS232', baud, timeout=tm) as paros:
                cmd = '*9900EW*9900BL=0\r\n'
                cmd = cmd.encode()
                paros.write(cmd)
                buf = paros.read(32).decode('utf-8')
                if '*9900EW' in buf and '*0001BL=0' in buf and '*9900BL=0' in buf:
                    break
                if baud == 300:
                    return False, 'No response from Paros sensor!'
        except Exception as e:
            return False, str(e)

    if baud == 9600:
        return True, 'Baud rate is already 9600!'
    try:
        with serial.Serial('RS232', baud, timeout=tm) as paros:
            cmd = '*9900EW*9900BR=9600\r\n'
            cmd = cmd.encode()
            paros.write(cmd)
            buf = paros.read(32).decode('utf-8')
            if '*9900EW' in buf and '9900BR=9600' in buf:
                return (
                    True,
                    'Baud rate changed from ' + str(baud) + ' to 9600!'
                )
    except Exception as e:
        return False, str(e)
    return (
        False,
        "Couldn't change Baud rate from " + str(baud) + ' to 9600!'
    )


def paros_get_data():
    """
    This function is used for retrieving paros data.
    """
    try:
        with serial.Serial('RS232', timeout=6) as paros:
            cmd = '*0100P3\r\n'
            cmd = cmd.encode()
            paros.write(cmd)
            buf = float(paros.readline().decode('utf-8')[5:].strip())
        return buf
    except Exception:
        return -99999.0


def paros_setup(paros_baud):
    """
    This function is used for setting up paros sensor.
    """
    if not paros_baud[0]:
        return False, [paros_baud[1]]
    response_list = [paros_baud[1]]
    err_count = 0
    command = 'TR=3840', 'UF=.6894757', 'PR=980', 'UN=0', 'MD=1', 'DP=4'
    try:
        with serial.Serial('RS232', timeout=1) as paros:
            for cmd in command:
                addr = '*0100'
                cmd = addr + 'EW' + addr + cmd + '\r\n'
                cmd = cmd.encode()
                paros.write(cmd)
                buf = paros.readline().decode('utf-8').strip()[1:]
                cmd = '0001' + cmd.decode('utf-8').strip()[12:]
                if cmd != buf:
                    if len(buf) == 0:
                        buf = 'no'
                    err_count += 1
                    response_list.append(
                        'Err: Got '
                        + buf + ' response when *0100'
                        + cmd[4:] + ' command was sent!'
                    )
                else:
                    response_list.append(buf + ' Successful!')
        if err_count == 0:
            return True, response_list
        return False, response_list
    except Exception as e:
        return False, [str(e)]


def ports_tag_message_append(flag, val, typ):
    """
    This function is used for appending ports tag messages.
    """
    msg = ''
    if typ == 1:
        msg += flag
        if -99999.0 not in val[:3]:
            msg += '{0}{1}{2}{3}{4}\r\n'.format('{:>11.3f}'.format(val[0]),
                                                '{:>9.3f}'.format(val[1]),
                                                '{:>10.0f}'.format(val[2]),
                                                '{:>10.1f}'.format(val[3]),
                                                '{:>10.1f}'.format(val[4]))
        else:
            msg += '  Data flagged as bad or missing\r\n'
    elif typ == 2:
        msg += flag
        if -99999.0 not in val:
            msg += '{0}{1}{2}\r\n'.format('{:>11.3f}'.format(val[0]),
                                          '{:>9.3f}'.format(val[1]),
                                          '{:>10.0f}'.format(val[2]))
        else:
            msg += '  Data flagged as bad or missing\r\n'
    elif typ == 3:
        msg += flag
        if -99999.0 not in val:
            msg += '{0}{1}{2}\r\n'.format('{:>11.1f}'.format(val[0]),
                                          '{:>9.0f}'.format(val[1]),
                                          '{:>10.1f}'.format(val[2]))
        else:
            msg += '  Data flagged as bad or missing\r\n'
    elif typ == 4 or typ == 6:
        msg += flag
        msg += '{:>11.1f}\r\n'.format(val) if val != -99999.0 else '  Data flagged as bad or missing\r\n'
    elif typ == 5:
        msg += flag
        msg += '{:>10.2f}\r\n'.format(val) if val != -99999.0 else ' Data flagged as bad or missing\r\n'
    elif typ == 7:
        msg += flag
        if -99999.0 not in val:
            msg += '{0}{1}\r\n'.format('{:>10.3f}'.format(val[0]),
                                       '{:>9.0f}'.format(val[1]))
        else:
            msg += '  Data flagged as bad or missing\r\n'
    elif typ == 8:
        msg += flag
        msg += '{:>10.3f}'.format(val) + '\r\n' if val != -99999.0 else ' data not available\r\n'
    elif typ == 9:
        for v in val:
            msg += flag
            msg += '{:>11.3f}'.format(v) + '\r\n' if v != -99999.0 else ' data not available\r\n'
    return msg


def ports_tag_message_formatter():
    """
    This function formats the data to a file for PORTS Tag transmission.
    """
    global goes_msg

    # Type 1
    # WL sensor with T1 and T2 = sns_goes, pri_temp, pri_rtemp,
    # pri_std_temp, pri_out_temp,
    # t1_temp, t2_temp, DPAS, GOES flag, redundant GOES Flag, sns_list,
    # format number, sns_label

    # Type 2
    # WL sensor sensor = sns_goes, pri_temp, pri_rtemp,
    # pri_std_temp, pri_out_temp,
    # DPAS, GOES flag, redundant GOES flag, sns_list,
    # format number, sns_label

    # Type 3
    # Triple met sensor (Wind) = sns_goes, sns1_temp, sns2_temp, sns3_temp,
    # DPAS, GOES flag, sns_list, format number, sns_label

    # Type 4
    # Single met sensor = sns_goes,
    # DPAS, GOES flag, format number, sns_label

    # Type 5
    # Two bytes flag single met sensor (Conductivity) = sns_goes,
    # DPAS, GOES flag, format number, sns_label

    # Type 6
    # Battery Voltage = sns_goes,
    # DPAS, GOES flag, format number, sns_label

    # Type 7
    # Two bytes flag double met sensor (Visibility) = sns_goes,
    # sns1_temp, sns2_temp,
    # DPAS, GOES flag, sns_list, format number, sns_label

    # Type 8
    # Data and sns offset = data_sns_goes,
    # format number, sns_label

    # Type 9
    # Tsunami water level = sns_goes,
    # DPAS, GOES flag, format number, sns_label

    all_sns_list = [[]] * 30
    if aq_active:
        # Type 1, GOES Flag: 1, Redundant GOES flag: >,
        # DPAS Code: A, Aquatrak Water Level
        all_sns_list[0] = [''] * 7 + ['A1', '1', '>', [], 1, all_sns[0]]
    else:
        # Type 2, GOES Flag: 1, Redundant GOES flag: >,
        # DPAS Code: A, Aquatrak Water Level
        all_sns_list[0] = [''] * 5 + ['A1', '1', '>', [], 2, all_sns[0]]

    # Type 2, GOES Flag: !, Redundant GOES flag: .,
    # DPAS Code: V, Shaft Angle Encoder 1
    all_sns_list[1] = [''] * 5 + ['V1', '!', '.', [], 2, all_sns[1]]

    # Type 2, GOES Flag: !, Redundant GOES flag: .,
    # DPAS Code: V, Shaft Angle Encoder 2
    all_sns_list[2] = [''] * 5 + ['V2', '!', '.', [], 2, all_sns[2]]

    # Type 2, GOES Flag: 8, Redundant GOES flag: #,
    # DPAS Code: Y, MWWL 1
    all_sns_list[3] = [''] * 5 + ['Y1', '8', '#', [], 2, all_sns[3]]

    # Type 2, GOES Flag: 8, Redundant GOES flag: #, DPAS Code: Y, MWWL 2
    all_sns_list[4] = [''] * 5 + ['Y2', '8', '#', [], 2, all_sns[4]]

    # Type 2, GOES Flag: %, Redundant GOES flag: ',
    # DPAS Code: N, Paroscientific Digiquartz #1 1
    all_sns_list[5] = [''] * 5 + ['N1', '%', "'", [], 2, all_sns[5]]

    # Type 2, GOES Flag: %, Redundant GOES flag: ',
    # DPAS Code: N, Paroscientific Digiquartz #1 2
    all_sns_list[6] = [''] * 5 + ['N2', '%', "'", [], 2, all_sns[6]]

    # Type 2, GOES Flag: %, Redundant GOES flag: ',
    # DPAS Code: T, Paroscientific Digiquartz #2 1
    all_sns_list[7] = [''] * 5 + ['T1', '&', '*', [], 2, all_sns[7]]

    # Type 2, GOES Flag: %, Redundant GOES flag: ',
    # DPAS Code: T, Paroscientific Digiquartz #2 2
    all_sns_list[8] = [''] * 5 + ['T2', '&', '*', [], 2, all_sns[8]]

    if ag_active:
        # Type 1, GOES Flag: (, Redundant GOES flag: ),
        # DPAS Code: Q, Airgap 1
        all_sns_list[9] = [''] * 7 + ['Q1', '(', ')', [], 1, all_sns[9]]
    else:
        # Type 2, GOES Flag: (, Redundant GOES flag: ),
        # DPAS Code: Q, Airgap 1
        all_sns_list[9] = [''] * 5 + ['Q1', '(', ')', [], 2, all_sns[9]]

    # Type 2, GOES Flag: (, Redundant GOES flag: ),
    # DPAS Code: Q, Airgap 2
    all_sns_list[10] = [''] * 5 + ['Q2', '(', ')', [], 2, all_sns[10]]

    # Type 3, GOES Flag: 3, DPAS Code: C, Wind 1
    all_sns_list[11] = [''] * 4 + ['C1', '3', [], 3, all_sns[11]]

    # Type 3, GOES Flag: 3, DPAS Code: C, Wind 2
    all_sns_list[12] = [''] * 4 + ['C2', '3', [], 3, all_sns[12]]

    # Type 4, GOES Flag: 4, DPAS Code: D, Air Temp 1
    all_sns_list[13] = ['', 'D1', '4', 4, all_sns[13]]

    # Type 4, GOES Flag: 4, DPAS Code: D, Air Temp 2
    all_sns_list[14] = ['', 'D2', '4', 4, all_sns[14]]

    # Type 4, GOES Flag: 5, DPAS Code: D, Water Temp 1
    all_sns_list[15] = ['', 'E1', '5', 4, all_sns[15]]

    # Type 4, GOES Flag: 5, DPAS Code: D, Water Temp 2
    all_sns_list[16] = ['', 'E2', '5', 4, all_sns[16]]

    # Type 4, GOES Flag: 6, DPAS Code: F, Barometric Pressure 1
    all_sns_list[17] = ['', 'F1', '6', 4, all_sns[17]]

    # Type 4, GOES Flag: 6, DPAS Code: F, Barometric Pressure 2
    all_sns_list[18] = ['', 'F2', '6', 4, all_sns[18]]

    # Type 5, GOES Flag: 6, DPAS Code: G, Conductivity 1
    all_sns_list[19] = ['', 'G1', '-7', 5, all_sns[19]]

    # Type 5, GOES Flag: 6, DPAS Code: G, Conductivity 2
    all_sns_list[20] = ['', 'G2', '-7', 5, all_sns[20]]

    # Type 4, GOES Flag: 9, DPAS Code: R, Relative Humidity 1
    all_sns_list[21] = ['', 'R1', '9', 4, all_sns[21]]

    # Type 4, GOES Flag: 9, DPAS Code: R, Relative Humidity 2
    all_sns_list[22] = ['', 'R2', '9', 4, all_sns[22]]

    # Type 7, GOES Flag: -O, DPAS Code: O, Visibility
    all_sns_list[23] = [''] * 3 + ['O1', '-O', [], 7, all_sns[23]]

    # Type 6, GOES Flag: =, DPAS Code: M, Aux Battery Voltage 1
    all_sns_list[24] = ['', 'M1', '=', 6, all_sns[24]]

    # Type 6, GOES Flag: <, DPAS Code: L, Battery Voltage 1
    all_sns_list[25] = ['', 'L1', '<', 6, all_sns[25]]

    # Type 2, GOES Flag: 2, Redundant GOES flag: ",
    # DPAS Code: B, Backup Water Level
    all_sns_list[26] = [''] * 5 + ['B1', '2', '"', [], 2, all_sns[26]]

    # Type 6, GOES Flag: <, DPAS Code: L, Battery Voltage 2
    all_sns_list[27] = ['', 'L1', '<', 6, all_sns[27]]

    # Type 8, SNS and DAT
    all_sns_list[28] = ['', 8, all_sns[28]]

    # Type 9, GOES Flag: T, DPAS Code: U, Tsunami
    all_sns_list[29] = ['', 'U1', 'T', 9, all_sns[29]]
    pri_year = '{:04d}'.format(add_sns[0].year)
    pri_month = '{:02d}'.format(add_sns[0].month)
    pri_day = '{:02d}'.format(add_sns[0].day)
    pri_date = '{0}/{1}/{2}'.format(pri_month, pri_day, pri_year)
    pri_hour = '{:02d}'.format(add_sns[0].hour)
    pri_minute = '{:02d}'.format(add_sns[0].minute)
    pri_second = '{:02d}'.format(add_sns[0].second)
    time_tag_goes = (
            '0' + add_sns[0].get_encoded_sutron_day()
            + add_sns[0].get_encoded_hour()
    )
    min_goes = add_sns[0].get_encoded_minute()
    sys_goes = '@@'
    ports_tag_msg = 'NOS {0} {1} {2}:{3}:{4}\r\n'.format(
        station_id,
        pri_date,
        pri_hour,
        pri_minute,
        pri_second
    )
    pri_goes = ''
    pri_sns_check = ''
    for asn in add_sns:
        for asl in range(len(all_sns_list)):
            if asn.label in all_sns_list[asl][-1]:

                # pri_goes is a placeholder for the primary sensor's goes message
                # pri_sns_check is a placeholder for the primary sensor's first label
                # ports_tag_msg holds the PORTS TAG's string as it's being built
                # all_sns_list[] is updating but its mutable so no need to return it
                # asn is for sensor class and is not updating
                (
                    pri_goes,
                    pri_sns_check,
                    ports_tag_msg
                ) = upd_meas(
                    pri_goes,
                    pri_sns_check,
                    ports_tag_msg,
                    all_sns_list[asl],
                    asn
                )

    # Remove SNS and DAT from all_sns_list and move its GOES message
    dat_sns_goes = all_sns_list.pop(-2)[0]

    goes_msg = 'P{0}{1}{2}{3}{4}{5}'.format(
        station_id,
        dat_sns_goes,
        sys_goes,
        min_goes,
        time_tag_goes,
        pri_goes
    )
    for asl in all_sns_list:
        goes_msg += asl[0]
    ports_tag_msg += '\r\nREPORT COMPLETE\r\n'
    with open('p3', 'w') as f:
        f.write(ports_tag_msg)


def pseudo_encoder(integer_value, byte_number, sign=False, decimal_places=0):
    """
    Pseudobinary encoder function converts decimal numbers to pseudobinary.
    :param integer_value: Decimal number.
    :param byte_number: Number of bytes.
    :param sign: Applies 2's compliment when sign is True.
    :param decimal_places: Use this shift.
    :return: Pseudobinary b format.
    """
    bits = 6 * byte_number
    integer_value = int(integer_value * 10 ** decimal_places)
    if sign:
        upper_limit = 2 ** bits // 2 - 1
        lower_limit = -(upper_limit + 1)
    else:
        upper_limit = 2 ** bits - 1
        lower_limit = 0
    integer_value = min(max(lower_limit, integer_value), upper_limit)
    if integer_value < 0:
        integer_value = 2 ** bits + integer_value
    binary_output = ''
    for by_num in range(byte_number):
        temp_value = integer_value >> 6 * by_num & 63
        temp_value += 64 if temp_value != 63 else 0
        binary_output = chr(temp_value) + binary_output
    return binary_output


def remove_old_samples(data_list, log_duration, end_time):
    """
    Keep only samples between [wv_sta, wv_end].
    Returns a NEW list (does not mutate one_hz_list).
    """
    n = len(data_list)
    sta_time = end_time - log_duration
    i = 0
    while i < n and data_list[i][0] < sta_time:
        i += 1
    j = n
    while j > i and data_list[j - 1][0] > end_time:
        j -= 1
    return data_list[i:j]


def return_vis_status():
    """
    This function returns visibility status.
    """
    vis_cmd('OPEN*').decode('utf-8', 'ignore') + '\r\n'
    v_status = vis_cmd('STATUS').decode('utf-8', 'ignore') + '\r\n'
    vis_cmd('CLOSE').decode('utf-8', 'ignore') + '\r\n'
    idx = v_status.find('PWD STATUS')
    if idx != -1:
        v_status = v_status[idx:]
    vis_msg = '{0}\r\n{1}\r\n'.format(sl3_datetime(), v_status)
    with open('vs', 'w') as f:
        f.write(vis_msg)
    return vis_msg


def sl3_datetime(date_time=None, tup=False):
    """
    This function changes the time and date as string or tuple
    -> ("YYYY/MM/DD HH:MM:SS").
    :return: Date and time.
    """
    if date_time is None:
        date_time = utime.localtime()[0:6]
    else:
        date_time = utime.localtime(date_time)[0:6]
    str_date = '{:04d}/{:02d}/{:02d}'.format(
        date_time[0],
        date_time[1],
        date_time[2]
    )
    str_time = '{:02d}:{:02d}:{:02d}'.format(
        date_time[3],
        date_time[4],
        date_time[5]
    )
    if tup:
        return str_date, str_time
    return str_date + ' ' + str_time


def sort_sns_list(list_item, list_label, old_list, new_list):
    """
    This function helps with sorting the sensor list in order.
    """
    for lst in list_item:
        try:
            new_list.append(old_list.pop(list_label.index(lst)))
            list_label.pop(list_label.index(lst))
        except ValueError:
            return


def status_message(msg):
    """
    This function updates the status file with status messages
    only when activated.
    """
    if gp5_status_msg >= 1:
        wait = utime.time()
        while utime.time() - wait < 10:
            if lock_write.acquire(False):
                sl3_date, sl3_time = sl3_datetime(tup=True)
                with open(
                        '/sd/status_log/status'
                        + sl3_date.replace('/', '.')
                        + '.txt', 'a'
                ) as f:
                    f.write(
                        '{0} {1}: {2}\r\n'.format(
                        sl3_date,
                        sl3_time, msg
                        )
                    )
                lock_write.release()
                print(msg)
                return
            utime.sleep_ms(50)

    else:
        print(msg)
        return


def sutron_day_calc(julian_day, year):
    """
    Sutron day counts the number of days from 12/31/1984 till date.
    Uses modulo 4096.
    :param julian_day: Julian day (between 1 and 365 or 366(ly))
    :param year: Year (yyyy)
    :return: Sutron day.
    """

    year -= 1985  # subtract 1984 + 1 from current year
    # Subtracting 1 because this year is not over
    leap_year = year // 4
    return (365 * year + leap_year + julian_day) % 4096


def upd_meas(pri_goes, pri_sns_check, ports_out, sns_para, asn):
    """
    This function updates measurements based on
    various formats and conditions.
    """
    if sns_para[-2] in (1, 2):
        # Example condition 1:
        # ['', '', '', '', '', '', '', 'Q1', '(', ')', [], 1, ['AG', 'AGSTD', 'AGOUT', 'AGT1', 'AGT2']]
        # Example condition 2:
        # ['', '', '', '', '', 'Y1', '8', '#', [], 2, ('MWWL', 'MWSTD', 'MWOUT')]
        if asn.label == sns_para[-1][0]:
            sns_para[1] = asn.get_encoded_data()
            sns_para[2] = asn.get_encoded_redundant_data()
            if asn.meas_number == 'M1':
                pri_sns_check = asn.label
        elif asn.label == sns_para[-1][1]:
            sns_para[3] = asn.get_encoded_data()
        elif asn.label == sns_para[-1][2]:
            sns_para[4] = asn.get_encoded_data()
        if sns_para[-2] == 1:
            if asn.label == sns_para[-1][3]:
                sns_para[5] = asn.get_encoded_data()
            elif asn.label == sns_para[-1][4]:
                sns_para[6] = asn.get_encoded_data()
        sns_para[-3].append(asn.value)
        if len(sns_para[-3]) == 3:
            if pri_sns_check == sns_para[-1][0]:
                pri_goes = (
                        sns_para[-5]
                        + sns_para[1]
                        + sns_para[3]
                        + sns_para[4]
                        + sns_para[-4]
                        + sns_para[2]
                )
            else:
                sns_para[0] = (
                        sns_para[-5]
                        + sns_para[1]
                        + sns_para[3]
                        + sns_para[4]
                        + sns_para[-4]
                        + sns_para[2]
                )
            if sns_para[-2] == 2:
                ports_out += ports_tag_message_append(
                    sns_para[-6] + ' ' + sns_para[-5],
                    sns_para[-3],
                    sns_para[-2]
                )
        elif len(sns_para[-3]) == 5:
            if pri_sns_check == sns_para[-1][0]:
                pri_goes = (
                        pri_goes[:-4]
                        + sns_para[5]
                        + sns_para[6]
                        + pri_goes[-4:]
                )
            else:
                sns_para[0] = (
                        sns_para[0][:-4]
                        + sns_para[5]
                        + sns_para[6]
                        + sns_para[0][-4:]
                )
            ports_out += ports_tag_message_append(
                sns_para[-6] + ' ' + sns_para[-5],
                sns_para[-3],
                sns_para[-2]
            )

    elif sns_para[-2] == 3:
        # Example condition 3:
        # ['', '', '', '', 'C1', '3', [], 3, ('WS', 'WD', 'WG')]
        if asn.label == sns_para[-1][0]:
            sns_para[1] = asn.get_encoded_data()
            if asn.meas_number == 'M1':
                pri_sns_check = asn.label
        elif asn.label == sns_para[-1][1]:
            sns_para[2] = asn.get_encoded_data()
        elif asn.label == sns_para[-1][2]:
            sns_para[3] = asn.get_encoded_data()
        sns_para[-3].append(asn.value)
        if len(sns_para[-3]) == 3:
            if pri_sns_check == sns_para[-1][0]:
                pri_goes = (
                        sns_para[-4]
                        + sns_para[1]
                        + sns_para[2]
                        + sns_para[3]
                )
            else:
                sns_para[0] = (
                        sns_para[-4]
                        + sns_para[1]
                        + sns_para[2]
                        + sns_para[3]
                )
            ports_out += ports_tag_message_append(
                sns_para[4] + ' ' + sns_para[5],
                sns_para[-3],
                sns_para[-2]
            )

    elif sns_para[-2] in (4, 5):
        # Example condition 4:
        # ['', 'D1', '4', 4, ('AT',)]
        # Example condition 5:
        # ['', 'G1', '-7', 5, ('COND',)]
        if asn.meas_number == 'M1':
            pri_goes = sns_para[2] + asn.get_encoded_data()
        else:
            sns_para[0] = sns_para[2] + asn.get_encoded_data()
        ports_out += ports_tag_message_append(
            sns_para[1] + ' ' + sns_para[2],
            asn.value,
            sns_para[-2]
        )

    elif sns_para[-2] == 6:
        # Example condition 6:
        # ['', 'L1', '<', 6, ('BAT',)]
        if sns_para[1] == 'L1':
            l1 = ' '
        else:
            l1 = ''
        sns_para[0] = sns_para[2] + asn.get_encoded_data() + l1
        ports_out += ports_tag_message_append(
            sns_para[1] + ' ' + sns_para[2],
            asn.value,
            sns_para[-2]
        )

    elif sns_para[-2] == 7:
        # Example condition 7:
        # ['', '', '', 'O1', '-O', [], 7, ('VIS3M', 'ALV')]
        if asn.label == sns_para[-1][0]:
            sns_para[1] = asn.get_encoded_data()
            if asn.meas_number == 'M1':
                pri_sns_check = asn.label
        elif asn.label == sns_para[-1][1]:
            sns_para[2] = asn.get_encoded_data()
        sns_para[-3].append(asn.value)
        if len(sns_para[-3]) == 2:
            if pri_sns_check == sns_para[-1][0]:
                pri_goes = (
                        sns_para[-4]
                        + sns_para[1]
                        + sns_para[2]
                )
            else:
                sns_para[0] = (
                        sns_para[-4]
                        + sns_para[1]
                        + sns_para[2]
                )
            ports_out += ports_tag_message_append(
                sns_para[3] + ' ' + sns_para[4],
                sns_para[-3],
                sns_para[-2]
            )

    elif sns_para[-2] == 8:
        # Example condition 8:
        # ['', 8, ('SNS', 'DAT')]
        if asn.label == sns_para[-1][0]:
            sns_para[0] = asn.get_encoded_data()
        elif asn.label == sns_para[-1][1]:
            sns_para[0] = asn.get_encoded_data() + sns_para[0]
        ports_out += ports_tag_message_append(
            asn.label,
            asn.value,
            sns_para[-2]
        )

    elif sns_para[-2] == 9:
        # Example condition 9:
        # ['', 'U1', 'T', 9, ('AQTWL', 'MWTWL')]
        if asn.label[:2] == pri_sns_check[:2]:
            tsu = asn.get_encoded_tsunami()
            sns_para[0] = (
                    sns_para[2]
                    + tsu[0] + tsu[1] + tsu[2]
                    + tsu[3] + tsu[4] + tsu[5]
                    + tsu[6] + tsu[7] + tsu[8]
            )
            val = (
                asn.value,
                asn.value2,
                asn.value3,
                asn.value4,
                asn.value5,
                asn.value6
            )
            ports_out += ports_tag_message_append(
                sns_para[1],
                val, sns_para[-2]
            )
    return pri_goes, pri_sns_check, ports_out


def update_dqap(sns_class, sns_name, meas_obj, deadline_off, log_duration):
    """
    This function checks for updated raw values and updates the DQAP or mean.
    """
    data_time = meas_as_reading(meas_obj).time
    deadline = utime.time() + deadline_off
    loop_success = False
    twl_mean = None
    while utime.time() < deadline:
        if sns_class.dqap_lock.acquire(False):
            try:
                if data_time <= sns_class.dqap_update:
                    dqap_list = remove_old_samples(
                        sns_class.dqap_list, log_duration, data_time
                    )
                    if sns_name.startswith('DQAP') or sns_name == 'PAROS':
                        dqap_list[:] = [sub[1] for sub in dqap_list if sub[1] > 0.0]
                        (
                            sns_class.mean,
                            sns_class.std,
                            sns_class.out,
                            sns_class.count
                        ) = calc_dqap(dqap_list)
                        sns_class.mean_update = utime.time()
                    elif sns_name == 'TWL':
                        dqap_list[:] = [sub[1] for sub in dqap_list if sub[1] > 0.0]
                        twl_mean = calc_mean(dqap_list)
                    elif sns_name == 'VIS':
                        dqap_list[:] = [sub[1] for sub in dqap_list]
                        (
                            sns_class.vis3,
                            sns_class.vis10,
                            sns_class.alv
                        ) = calc_vis(dqap_list)
                        sns_class.mean_update = utime.time()

                    loop_success = True
                    break
            finally:
                sns_class.dqap_lock.release()
        utime.sleep_ms(50)
    if loop_success:
        if sns_name.startswith('DQAP') or sns_name == 'PAROS':
            if sns_class.mean == -99999.0:
                status_message('{0} mean unsuccessful!'.format(sns_name))
            else:
                status_message('{0} mean update successful!'.format(sns_name))
            return
        if sns_name == 'TWL':
            if twl_mean == -99999.0:
                status_message(
                    'Tsunami{0} mean update unsuccessful...'.format(
                        sns_class.object_name[-1]
                    )
                )
            else:
                status_message(
                    'Tsunami{0} mean update successful!'.format(
                        sns_class.object_name[-1]
                    )
                )
            return twl_mean
        if sns_name == 'VIS':
            if sns_class.vis3 == -99999.0:
                status_message('VIS3M update unsuccessful!')
            else:
                status_message('VIS3M update successful!')
            return
    if sns_name.startswith('DQAP') or sns_name == 'PAROS':
        status_message(
            '{0} mean update unsuccessful, waited too long...'.format(sns_name)
        )
        sns_class.reset()
        return
    if sns_name == 'TWL':
        status_message(
            'Tsunami{0} mean update unsuccessful, waited too long...'.format(
                sns_class.object_name[-1]
            )
        )
        return -99999.0
    if sns_name == 'VIS':
        status_message('VIS3M update unsuccessful, waited too long...')
        sns_class.reset()
        return


def val_tag(val=-99999.0):
    """
    This function tags the measurement as bad when the value is -99999.0.
    """
    if val == -99999.0:
        meas_make_invalid()
    return val


def vis_bitostr(bi_data):
    """
    This function converts integers representing ASCII values to string.
    """
    bi_data = [bd if 127 > bd > 31 else "nan" for bd in bi_data]
    dec_data = ''
    for bd in bi_data:
        if bd not in ("nan", 62):
            dec_data += chr(bd)
    return dec_data.strip()


def vis_cmd(command, vis_ini=True):
    """
    This function sends the command to the visibility sensor.
    """
    with serial.Serial('RS232') as vis:
        if vis_ini or command == 'STATUS':
            vis.timeout = 2
            vis.inter_byte_timeout = 1
        else:
            vis.timeout = 0.5
            vis.inter_byte_timeout = 0.1
        vis.reset_input_buffer()
        vis.reset_output_buffer()
        vis.write(bytes(command + '\r', 'utf-8'))
        response = vis.read(1024)
    return response


def vis_get_data():
    """
    This function gets the data from the visibility sensor.
    """
    vis_cr = chr(13)
    vis_enq = chr(5)
    vis_id_space = '  '
    vis_msg = '0'
    vis_data_cmd = vis_cr + vis_enq + 'PW' + vis_id_space + vis_msg
    vis = vis_cmd(vis_data_cmd, False)
    vis_out = [-99999.0] * 3
    try:
        vis_out[0] = int(vis[7:9])
    except ValueError:
        pass
    try:
        vis_out[1] = int(vis[9:16]) * 1e-3
    except ValueError:
        pass
    try:
        vis_out[2] = int(vis[16:22]) * 1e-3
    except ValueError:
        pass
    return vis_out


def vis_initialize():
    """
    This function initializes the visibility sensor.
    """
    vis_response = []
    cmd = ('OPEN', 'AMES 0 0', 'BAUD', 'CLOSE')
    vis_check = ('PWD', 'AUTOMATIC', 'BAUD', 'LINE')
    vis_success = True
    vis_cmd('CLOSE')
    for c, v in zip(cmd, vis_check):
        response = vis_bitostr(vis_cmd(c))
        if v in response:
            vis_response.append(response + ', Successful!')
        else:
            vis_success = False
            vis_response.append(c + ' command failed...')
    return vis_success, vis_response


@TASK
def init():
    """
    This task initializes the active sensors.
    """
    global add_sns, ag, ag_active, aq, aq_active, cnt_meas
    global curr_time, initialize, m1_sns_check, no_of_meas
    global paros_class, paros_init, vis_class, vis_init
    global dqap1_class, dqap1_init, dqap2_class, dqap2_init
    add_sns = []
    temp_sns = []
    temp_label = []
    cnt_meas = 0

    try:
        clear_dir()
        with open('p3', 'w') as f:
            f.write(
                '{0}: Initializing, please check again in a few minutes...\r\n'.format(
                    sl3_datetime()
                )
            )

        while no_of_meas < 65:
            response = setup_read('!M{0}'.format(no_of_meas))
            if 'Unknown command' in response:
                break
            no_of_meas += 1



        command_line('!file mkdir /sd/status_log/\r')
        status_message('Initializing configuration...')
        curr_time = int(utime.time())

        for sn in range(1, no_of_meas):
            if setup_read('!M{0} Active'.format(sn)) == 'On':
                sns_label = setup_read('!M{0} Label'.format(sn))
                if sn == 1:
                    m1_sns_check = sns_label

                # Check for DQAP measurement
                if sns_label.startswith('DQAP'):
                    dqap_temp = sns_label
                    sns_label = 'DQAP'
                    dqap_name = dqap_temp[:5]
                    if dqap_name in ('DQAP1', 'DQAP2'):
                        if len(dqap_temp) > 6:
                            try:
                                dqap_temp = int(
                                    float(dqap_temp.split('_')[1])
                                )
                                if dqap_name == 'DQAP1':
                                    dqap1_class = DQAP(dqap_temp, 'DQAP1')
                                    dqap1_init = True

                                elif dqap_name == 'DQAP2':
                                    dqap2_class = DQAP(dqap_temp, 'DQAP2')
                                    dqap2_init = True
                                status_message(
                                    '{0} initialization successful!'.format(
                                        dqap_name
                                    )
                                )
                            except Exception:
                                status_message(
                                    '{0} failed to initialize. Please check measurement label...'.format(
                                        dqap_name
                                    )
                                )

                # Check for Visibility measurement
                if sns_label == 'VISRAW':
                    vis_class = VisData()
                    vis_init = True
                    init_vis()

                # Check for Paros measurement
                if sns_label == 'PAROSRAW':
                    paros_class = DQAP(180, 'PAROS')
                    paros_init = True
                    init_paros()

                # Check for other measurements
                if sns_label not in discard_sns:
                    cnt_meas += 1
                    if cnt_meas == 1 or sns_label in primary_sns:
                        temp_sns.append(PrimarySensor('M' + str(sn)))
                    elif sns_label in tsunami_sns:
                        temp_sns.append(TsunamiData('M' + str(sn)))
                    else:
                        temp_sns.append(SecondarySensor('M' + str(sn)))

        for ts in temp_sns:
            temp_label.append(ts.label)

        # Check for T1 and T2
        if 'AQT1' in temp_label and 'AQT2' in temp_label:
            aq += aq_temp[1]
            aq_active = True
        else:
            if 'AQT1' in temp_label or 'AQT2' in temp_label:
                cnt_meas -= 1
        if 'AGT1' in temp_label and 'AGT2' in temp_label:
            ag += ag_temp[1]
            ag_active = True
        else:
            if 'AGT1' in temp_label or 'AGT2' in temp_label:
                cnt_meas -= 1

        # Check for errors in measurement labels. Will not initialize if not set correctly
        if check_label_err(temp_label):
            return

        # This loop appends primary sensor to the add_sns list
        one_meas = True
        for psg in primary_sns_grp:
            if psg[0] == temp_label[0]:
                one_meas = False
                sort_sns_list(psg, temp_label, temp_sns, add_sns)
                break
        if one_meas:
            add_sns.append(temp_sns.pop(0))
            temp_label.pop(0)

        # This loop appends the rest of the active sensors to the add_sns list
        for als in all_sns[:28] + all_sns_remaining:
            if als[0] in temp_label:
                sort_sns_list(als, temp_label, temp_sns, add_sns)
        ports_tag_message_formatter()
        initialize = True
        status_message('Initialization complete!')
    except Exception as e:
        status_message('Could not initialize...')
        status_message('Error: ' + str(e))


@TASK
def update_data():
    """
    This task updates the sensor objects.
    """
    if not initialize:
        status_message('System has not initialized yet!')
        return
    global curr_time
    status_message('Updating all data...')

    curr_time = int(utime.time())

    for i in range(cnt_meas):
        if i == 0 or add_sns[i].label in primary_sns:
            add_sns[i].update_primary_data()
        else:
            if add_sns[i].label in tsunami_sns:
                add_sns[i].update_tsunami_data()
            else:
                add_sns[i].update_secondary_data()
    ports_tag_message_formatter()
    status_message('Updates successful!')


@TASK
def vis_status():
    """
    This task updates visibility status and saves in a file named vs.
    """
    vis_class.return_vis_status()


@TXFORMAT
def goes_message(_):
    """
    This transmission function returns the GOES message for transmission.
    """
    if not initialize:
        return 'System has not initialized yet!'
    status_message('Transmitting GOES message...')
    good_goes_message = goes_msg
    print("1", good_goes_message)
    tx_battery = max(
        float(command_line('!BATT DURING LOAD\r').strip()), 9.5
    )
    tx_battery = round((tx_battery - 9.5) * 10)
    tx_battery = pseudo_encoder(tx_battery, 1)
    idx = good_goes_message.rfind(' ')
    cnt = good_goes_message.count(' ')
    if idx == -1:
        pass
    else:
        good_goes_message = (
                good_goes_message[:idx + 1]
                + tx_battery
                + good_goes_message[idx + 1:]
        )
        print("2", good_goes_message)
        if cnt == 2:
            good_goes_message = good_goes_message.replace(' ', '', 1)
    status_message('GOES transmission successful!')
    print("3", good_goes_message)
    return good_goes_message


@MEASUREMENT
def dqap1_raw(data_recent):
    """
    This measurement builds the data list.
    """
    meas_obj = data_recent
    data_recent = live_data()
    if not dqap1_init:
        return val_tag()
    dqap1_class.build_list(data_recent, meas_obj)
    return val_tag(data_recent)


@MEASUREMENT
def dqap1_mean(meas_obj):
    """
    This measurement function returns DQAP1 mean.
    """
    if not dqap1_init:
        status_message(
            'DQAP1 mean measurement has not initialized yet...'
        )
        return val_tag()
    return val_tag(dqap1_class.update_dqap(meas_obj))


@MEASUREMENT
def tsunami1_mean(meas_obj):
    """
    This measurement function returns mean from the
    last 60 seconds samples.
    """
    if not dqap1_init:
        status_message(
            'Tsunami1 mean measurement has not initialized yet...'
        )
        return val_tag()
    return val_tag(dqap1_class.tsunami_avg(meas_obj))


@MEASUREMENT
def dqap1_std(_):
    """
    This measurement function returns DQAP1 std.
    """
    if not dqap1_init:
        return val_tag()
    return val_tag(dqap1_class.dqap_parameter(0))


@MEASUREMENT
def dqap1_out(_):
    """
    This measurement function returns DQAP1 outliers.
    """
    if not dqap1_init:
        return val_tag()
    return val_tag(dqap1_class.dqap_parameter(1))


@MEASUREMENT
def dqap1_cnt(_):
    """
    This measurement function returns DQAP1 count.
    """
    if not dqap1_init:
        return val_tag()
    return val_tag(dqap1_class.dqap_parameter(2))


@MEASUREMENT
def dqap2_raw(data_recent):
    """
    This measurement builds the DQAP2 list.
    """
    meas_obj = data_recent
    data_recent = live_data()
    if not dqap2_init:
        return val_tag()
    dqap2_class.build_list(data_recent, meas_obj)
    return val_tag(data_recent)


@MEASUREMENT
def dqap2_mean(meas_obj):
    """
    This measurement function returns DQAP2 mean.
    """
    if not dqap2_init:
        status_message(
            'DQAP2 mean measurement has not initialized yet...'
        )
        return val_tag()
    return val_tag(dqap2_class.update_dqap(meas_obj))


@MEASUREMENT
def tsunami2_mean(meas_obj):
    """
    This measurement function returns mean
    from the last 60 seconds samples.
    """
    if not dqap2_init:
        status_message(
            'Tsunami2 mean measurement has not initialized yet...'
        )
        return val_tag()
    return val_tag(dqap2_class.tsunami_avg(meas_obj))


@MEASUREMENT
def dqap2_std(_):
    """
    This measurement function returns DQAP2 std.
    """
    if not dqap2_init:
        return val_tag()
    return val_tag(dqap2_class.dqap_parameter(0))


@MEASUREMENT
def dqap2_out(_):
    """
    This measurement function returns DQAP2 outliers.
    """
    if not dqap2_init:
        return val_tag()
    return val_tag(dqap2_class.dqap_parameter(1))


@MEASUREMENT
def dqap2_cnt(_):
    """
    This measurement function returns DQAP2 count.
    """
    if not dqap2_init:
        return val_tag()
    return val_tag(dqap2_class.dqap_parameter(2))


@MEASUREMENT
def paros_raw(meas_obj):
    """
    This measurement collects the raw paros samples from the sensor.
    """
    if not paros_init:
        return val_tag()
    # data_recent = paros_get_data()
    data_recent = live_data()
    paros_class.build_list(data_recent, meas_obj)
    return val_tag(data_recent)


@MEASUREMENT
def paros_mean(meas_obj):
    """
    This measurement function returns Paros mean.
    """
    if not paros_init:
        return val_tag()
    return val_tag(paros_class.update_dqap(meas_obj))


@MEASUREMENT
def paros_std(_):
    """
    This measurement function returns Paros standard deviation.
    """
    if not paros_init:
        return val_tag()
    return val_tag(paros_class.dqap_parameter(0))


@MEASUREMENT
def paros_out(_):
    """
    This measurement function returns the Paros out.
    """
    if not paros_init:
        return val_tag()
    return val_tag(paros_class.dqap_parameter(1))


@MEASUREMENT
def paros_cnt(_):
    """
    This measurement function returns the Paros counts.
    """
    if not paros_init:
        return val_tag()
    return val_tag(paros_class.dqap_parameter(2))


@MEASUREMENT
def vis_raw(meas_obj):
    """
    This measurement collects the raw visibility samples from the sensor.
    """
    if not vis_init:
        return val_tag()
    # data_recent = vis_get_data()
    data_recent = round(uniform(0,1)), live_data(), live_data()
    vis_class.build_list(data_recent, meas_obj)
    return val_tag(data_recent[1])


@MEASUREMENT
def vis_3m(meas_obj):
    """
    This measurement returns VIS3M.
    """
    if not vis_init:
        return val_tag()
    return val_tag(vis_class.update_vis(meas_obj))


@MEASUREMENT
def vis_10m(_):
    """
    This measurement returns VIS10M.
    """
    if not vis_init:
        return val_tag()
    return val_tag(vis_class.vis_parameter(0))


@MEASUREMENT
def alv(_):
    """
    This measurement returns ALV.
    """
    if not vis_init:
        return val_tag()
    return val_tag(vis_class.vis_parameter(1))
