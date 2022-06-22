#!/usr/bin/python
# if you get an ImportError, type this on the command line:
import os
import datetime
# DO: pip install pandas
import pandas as pd

def read_file(root, filename):
        lines = []
        full_path = root + '/' + filename
        with open(full_path, 'r') as file:
                for line in file:
                        if line.startswith('!'):
                                continue
                        else:
                                lines.append(line.split())
        return lines


def convert_to_dataframe(data1):
        new_df = pd.DataFrame(data=data1)
        return new_df


def write_to_csv(df):
        current_time = datetime.datetime.now()
        current_time = current_time.strftime("%H_%M_%S")
        filename = 'spreadsheet_formatted_{}'.format(current_time) + ".csv"
        df.to_csv(filename, header=False, index=False)


root = os.getcwd()
filename = 'crystal.txt'

data = read_file(root, filename)
df = convert_to_dataframe(data)
