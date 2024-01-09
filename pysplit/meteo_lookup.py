from ftplib import FTP
from datetime import datetime
from math import ceil, floor
from pathlib import Path
from calendar import month_name
from dateutil.relativedelta import *

""" meteo_handler provides functionality for interfacing with NOAA ARL datasets

For descriptions of the datasets and the protocol for downloading from the ftp server:
https://www.ready.noaa.gov/archives.php

Additional information on the datasets is given by:
https://www.ready.noaa.gov/documents/ppts/Cheat_Sheet_2020.pdf
"""

# arl_datasets only includes datasets that are actively updated
arl_datasets = {
    'nams': {
        'start': datetime(2010, 1, 1),
        'file_name': "%Y%m%d_hysplit.t00z.namsa",
        'folder': 'nams/',
        'time_step': relativedelta(days=+1),
        'file_size': 1200},
    'nams-AK': {
        'start': datetime(2010, 1, 1),
        'file_name': "%Y%m%d_hysplit.t00z.namsa.AK",
        'folder': 'nams/',
        'time_step': relativedelta(days=+1),
        'file_size': 712},
    'nams-HI': {
        'start': datetime(2010, 1, 1),
        'file_name': "%Y%m%d_hysplit.t00z.namsa.HI",
        'folder': 'nams/',
        'time_step': relativedelta(days=+1),
        'file_size': 481},
    'gdas': {
        'start': datetime(2004, 12, 1),
        'file_name': "gdas1.%x%Y.w%w",
        'folder': 'gdas/',
        'time_step': relativedelta(days=+1),
        'file_size': 571,
        'x': lambda x: month_name[x.month][0:3]},
    'gfs0p25': {
        'start': datetime(2004, 12, 1),
        'file_name': "%Y%m%d_gfs0p25",
        'folder': 'gfs0p25/',
        'time_step': relativedelta(days=+1),
        'file_size': 2700},
    'reanalysis': {
        'start': datetime(1948, 1, 1),
        'file_name': "RP%Y%m.gbl",
        'folder': 'reanalysis/',
        'time_step': relativedelta(months=+1),
        'file_size': 117},
    'hrrr': {
        'start': datetime(2019, 6, 1),
        'file_name': "%Y%m%d_%x_hrrr",
        'folder': 'hrrr/',
        'time_step': relativedelta(hours=+6),
        'file_size': 3200,
        'x': lambda x: ['00-05', '06-11', '12-17', '18-23'][floor(x.hour / 6)]},
    'nam12': {
        'start': datetime(2007, 5, 1),
        'file_name': "%Y%m%d_nam12",
        'folder': 'nam12/',
        'time_step': relativedelta(days=+1),
        'file_size': 275},
    'wrf27km-avg': {
        'start': datetime(1980, 1, 1),
        'file_name': 'wrfout_d01_%Y%m%d.ARL',
        'folder': 'wrf27km/avg/%Y/',
        'time_step': relativedelta(days=+1),
        'file_size': 275},
    'wrf27km-inst': {
        'start': datetime(1980, 1, 1),
        'file_name': 'wrfout_d01_%Y%m%d.ARL',
        'folder': 'wrf27km/inst/%Y/',
        'time_step': relativedelta(days=+1),
        'file_size': 210},
}


def week_of_month(date):
    """ Returns week of the month from a datetime variable.
    This is currently only required for the gdas dataset.

    Parameters
    date - datetime variable

    Returns
    Week of the month, ranges from 1 to 6

    Examples
    datetime(2023, 12, 31) → 6
    datetime(2024, 01, 01) → 1
    """
    first_day = datetime(date.year, date.month, 1)
    day_of_week = [x for x in range(1, 7)]
    day_of_week.append(0)
    offset = day_of_week[first_day.weekday()]
    return ceil((date.day + offset) / 7)


def parse_format(input_string):
    """ Helper function to find each delimiter % that exists in a string

    Parameters
    ----------
    input_string = string which may include % as a format code.

    Returns
    -------
    list of unique format codes. These are the first character that follows a %

    Examples
    -------
    parse_format('%f%m%y%Y%f') → ['f', 'm', 'Y']
    """
    format_code = []
    while True:
        idx = input_string.find('%')
        if idx == -1:
            break
        format_code.append(input_string[idx:(idx + 2)])
        input_string = input_string[idx + 2:]
    return list(set(format_code))


def format_met_string(date, file_name, **kwargs):
    """ Function to format strings using standard datetime delimiters

    Parameters
    date - datetime value
    file_name - string with delimiters marked by '%' which are to be
                replaced by values from the provided date.
    kwargs - specify delimiter and how it is formatted

    Returns
    file_name with inserted date values

    Examples
    d = datetime(2020, 1, 1)
    fname = 'text%y%m%d.cfg'
    format_met_string(d, fname)
        → 'text20200101.cfg'

    import calendar
    fname = 'other_%x.mon'
    func = lambda x: month_name[x.month]
    format_met_string(d, fname, x=func)
        → 'other_January.mon'
    """
    formatter = {
        '%Y': lambda x: str(x.year),
        '%y': lambda x: str(x.year)[2:],
        '%m': lambda x: str(x.month).zfill(2),
        '%w': lambda x: str(week_of_month(x)),
        '%d': lambda x: str(x.day).zfill(2),
        '%h': lambda x: str(x.hour).zfill(2),
    }
    for code in kwargs:
        formatter[f'%{code}'] = kwargs[code]

    codes = parse_format(file_name)
    for code in codes:
        file_name = file_name.replace(code, formatter[code](date))

    return file_name


def arl_file(date, dataset):
    """ Finds a single meteo file from the specified dataset that contains met
    data for the date requested.

    Parameters
    ----------
    date - datetime object
    dataset - string that is one of the ARL datasets. Must be a valid key in the
     arl_datasets dictionary.

    Returns
    -------
    string that is a single meteo data file name from the ARL

    Examples
    -------
    d = datetime(2020, 1, 1)
    arl_file(d, 'gdas') → 'gdas1.Jan2020.w1'
    """
    name_format = arl_datasets[dataset]['file_name']
    try:
        x = arl_datasets[dataset]['x']
        file_name = format_met_string(date, name_format, x=x)
    except KeyError:
        file_name = format_met_string(date, name_format)

    return file_name


def arl_span(dates, dataset):
    """ Finds the met data folder/filename for a given date in a given dataset within the NOAA ARL ftp server

    Parameters
    date_range - list of datetime variable that define the start and end
                 of the required time span. Ordering doesn't matter.
                 Also accepts a single datetime value. This should result
                 in a single returned filename.
    dataset - name of the specified ARL dataset. Valid dataset names
              are defined above under _valid_datasets.

    Returns
    list of all filenames for the given range using ARL convention for the specified dataset.
    This also returns the folder location for each met file.

    Examples
    date_range = [datetime(2020, 1, 1), datetime(2020, 1, 4)]
    meteo_name(date_range, 'nams') →
        ['20200101_hysplit.t00z.namsa',
         '20200102_hysplit.t00z.namsa',
         '20200103_hysplit.t00z.namsa',
         '20200104_hysplit.t00z.namsa']
    """
    if isinstance(dates, list):
        if dates.__len__() != 2:
            raise ValueError('Input date range is not two dates.')
        else:
            for date in dates:
                if not isinstance(date, datetime):
                    raise ValueError('Input must be datetime.')
        start_date = min(dates)
        end_date = max(dates)
    else:
        raise ValueError('Invalid date input. Must be list of datetime.')

    meteo_names = [arl_file(start_date, dataset)]
    while start_date < end_date:
        start_date += arl_datasets[dataset]['time_step']
        next_file = arl_file(start_date, dataset)
        if next_file != meteo_names[-1]:
            meteo_names.append(next_file)

    return meteo_names


def arl_list(dates, dataset):
    """ Similar to arl_span, but only finds files for a specific list.
    This is useful for when finding which files to download from the
    noaa arl ftp server without downloading excessive amounts of data.

    Parameters
    ----------
    dates - list of datetime variables
    dataset - name of the specified ARL dataset. Valid dataset names
              are defined above under _valid_datasets.

    Returns
    -------
    list of all filenames for the given range using ARL convention for the specified dataset.
    This also returns the folder location for each met file.

    Examples
    -------
    date_range = [datetime(2020, 1, 1), datetime(2020, 1, 4)]
    meteo_name(date_range, 'nams') →
        ['20200101_hysplit.t00z.namsa',
         '20200104_hysplit.t00z.namsa']
    """
    meteo_names = [arl_file(dates[0])]
    for date in dates[1:]:
        next_file = arl_file(date, dataset)
        if next_file != meteo_names[-1]:
            meteo_names.append(next_file)

    return meteo_names


def download_meteo_data(save_folder, meteo_files, dataset):
    """ Bulk download of met data from NOAA ARL ftp server.

    Parameters
    save_folder - Folder used on the local machine where met data is saved.
    date_range - datetime variable or list of datetime variables. The function
                 will automatically pick out the earliest and latest values and
                 compile the entire range of available data.
    dataset - name of the specified ARL dataset.

    Returns
    downloads all the met data files requested by the user.

    Warnings
    Met data files are quite large. Using this function without care could
    consume a large chunk of disk space if used carelessly.

    Examples
    save_folder = 'd:/met data/'
    start_date = datetime(1979, 1, 1)
    end_date = datetime(1979, 12, 1)
    date_range = [start_date, end_date]
    download_meteo_data(save_folder, date_range, 'reanalysis')
    → Files are downloaded directly to the specified directory
    """
    # look up expected download size
    size = arl_datasets[dataset]['file_size'] * meteo_files.__len__()
    if size > 1000:
        total = str(round(size / 1024, 2)) + 'GB'
    else:
        total = str(size) + 'MB'
    if input(f"Requested files will take up {total} of disk space. Do You Want To Continue? [y/n]: ") != "y":
        return
    site = 'ftp.arl.noaa.gov'
    user = 'anonymous'
    meteo_folder = arl_datasets[dataset]['folder']
    arl_ftp = FTP(host=site, user=user)
    arl_ftp.cwd(f'/archives/{meteo_folder}')
    for meteo_file in meteo_files:
        path = Path(save_folder + meteo_folder + meteo_file)
        if not path.is_file():
            print(f'downloading {meteo_file}')
            with open(save_folder + meteo_folder + meteo_file, 'wb') as file:
                arl_ftp.retrbinary(f'RETR {meteo_file}', file.write)

    arl_ftp.close()


if __name__ == '__main__':
    # Test arl_file database
    now = datetime(2024, 1, 1)
    for i in arl_datasets.keys():
        arl_file(now, i)
    print('arl_file passed.')

    # Test arl_span
    future = now + relativedelta(months=+12)
    date_range = [now, future]
    len = []
    for i in arl_datasets.keys():
        len.append(arl_span(date_range, i).__len__())
    print('arl_span passed.')
    print(len)

    # Test arl_list
    start_date = datetime(1978, 1, 1)
    end_date = datetime(1987, 5, 31)
    date_range = [start_date, end_date]
    reanalysis_files = arl_span(date_range, 'reanalysis')
    download_meteo_data('D:/met data/', reanalysis_files, 'reanalysis')
