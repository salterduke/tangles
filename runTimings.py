import os

for i in range(5):
    os.system("cp const_150.txt const.txt")
    print("Before VY 150, run: {}".format(i))
    os.system("python main.py VY")
    print("Before YWS 150, run: {}".format(i))
    os.system("python main.py VY")
    os.system("cp const_200.txt const.txt")
    print("Before VY 200, run: {}".format(i))
    os.system("python main.py VY")
    print("Before YWS 150, run: {}".format(i))
    os.system("python main.py VY")