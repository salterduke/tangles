import os

for i in range(5):
    os.system("cp const_all.txt const.txt")
    print("Before VY all, run: {}".format(i))
    os.system("python main.py VY")
    print("Before YWS all, run: {}".format(i))
    os.system("python main.py YWS")

# for i in range(3):
#     os.system("cp const_150.txt const.txt")
#     print("Before VY 150, run: {}".format(i))
#     os.system("python main.py VY")
#     print("Before YWS 150, run: {}".format(i))
#     os.system("python main.py YWS")
#     os.system("cp const_200.txt const.txt")
#     print("Before VY 200, run: {}".format(i))
#     os.system("python main.py VY")
#     print("Before YWS 200, run: {}".format(i))
#     os.system("python main.py YWS")