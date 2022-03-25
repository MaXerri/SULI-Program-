import os
import pandas as pd
import xlwt
from tempfile import TemporaryFile

def read_file(root, filename):
    lines = []
    full_path = root + '/' + filename
    with open(full_path, 'r') as file:
        for line in file:
            lines.append(line.split())
    return lines


def write_to_excel(data1):
    book = xlwt.Workbook()
    sheets = ['sheet1']

    for sheet, data in zip(sheets, [data1]):
        ws = book.add_sheet(sheet)

        for row, row_value in enumerate(data):
            for col, col_value in enumerate(row_value):
                ws.write(row, col, col_value)

    name = 'output.xls'
    book.save(name)


root = os.getcwd()
filename = 'sampledataset.txt'

data1 = read_file(root, filename)
write_to_excel(data1)
