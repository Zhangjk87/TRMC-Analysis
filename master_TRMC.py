import os
import configparser


def main():

    config = configparser.RawConfigParser()
    config.read('config.ini')
    filename = config['Constants']['filename']
    os.system("respeak3.py {0}".format(filename))
    os.system("base2.py {0}".format(filename))

if __name__ == '__main__':
    main()
