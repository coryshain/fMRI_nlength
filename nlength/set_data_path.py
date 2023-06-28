import os

data_path = input('Enter path to data repository (available at https://osf.io/fduve/) >>> ')
with open('data_path.txt', 'w') as f:
    f.write(data_path)