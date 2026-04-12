import csv
import sys
import numpy as np
from geo_helper import GeoHelper as gh
import matplotlib.pyplot as plt

def main():
    labels = [] 
    RT = []

    source = input("Enter source: 1 - Wells Simpson, 2 - large set, 3 - quadratic test")
    if source == "1":
        print("parameters from csv from Wells Simpson paper (edited)")
        file_path = "data/WellsSimpson2001.csv"
        title = file_path;
        print(f"Selected file: {file_path}")
        with open(file_path, mode='r', newline='', encoding='utf-8') as file:
            reader = csv.DictReader(file)
            for row in reader:
                if row["R"] != "" and row["T"] != "" and row["Use"] == "Y":
                    RT.append([float(row["R"]),float(row["T"])])
                    labels.append(row["Tag"])

    if source == "2":
        print("test parameters from csv using published paleo")
        file_path = "data/Unit,Sub-unitSite,Age_Ma,Rotation_d.csv"
        title = file_path;
        print(f"Selected file: {file_path}")
        with open(file_path, mode='r', newline='', encoding='utf-8') as file:
            reader = csv.DictReader(file)
            for row in reader:
                if row["Rotation_deg"] != "" and row["Age_Ma"] != "":
                    RT.append([float(row["Rotation_deg"]), float(row["Age_Ma"])])

    if source == "3":
        title = "quadratic fit test"
        print("test parameters from quadratic fit to graph:")
        a = -0.0088 # from quadratic fit
        b = 0.4766  # from quadratic fit
        w_0 = 1.19  # from Wells - Simpson Fig3B
        print("a = ", a)
        print("b = ", b)
        print("w_o = ", w_0)
        RT = createTestDataSet(a, b, w_0, 50, 10)

    #print ("RT: ", RT)
    c = calculateRateCoeff(RT, normalize=False)
    print("coeffeciante c: ", c)
    Erms = calculateRegressionRmsError(c, RT)
    print("Calibration residual rms: ", Erms)

    # Create the plots 
    # Paleomag data plot
    fig, ax1 = plt.subplots()
    Ri = col(RT, 0)
    Ti = col(RT, 1)
    ax1.scatter(Ti, Ri, label='Paleomag Data', color='blue')
    ax1.set_ylim(0)
    for i, txt in enumerate(labels):
        pad = 1
        ax1.annotate(txt, (Ti[i] + pad, Ri[i] + pad))
    ax1.set_ylabel('Total Rotation (Deg)', color='tab:blue')
    ax1.set_xlabel('T (Ma)')
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    # quadratic coefficients speed plot
    ax2 = ax1.twinx()
    T = np.linspace(0, 50, 11)
    w = calculateRotationSpeedFromCoeff(c, T)
    ax2.plot(T, w, color='orange', linestyle='--')
    ax2.set_ylim(0)
    ax2.set_xlabel('T (Ma)')
    ax2.set_ylabel('Quadratic Rot Rate (Deg / Ma)', color='tab:orange')

    # quadratic recalculation of total rotation
    Rc = calculateTotalRotationFromCoeff(c, T)
    ax1.plot(T, Rc, color='blue', linestyle='-')
    
    # 3. Add details
    plt.title("Gain Curve From " + title) 
    plt.grid(True)

    # 4. Display or save
    plt.show()  # Opens a window to show the plot
    # plt.savefig('my_plot.png') # Saves to file
    
    return

def col(maRTix, colNo):
    return [row[colNo] for row in maRTix]

def calculateRateCoeff(RTArray, normalize = False):
    T1 = T2 = T3 = T4 = T5 = T6 = RT1 = RT2 = RT3 = 0.0

    for RT in RTArray:
        T = RT[1]
        R = RT[0]
        # print (RT)
        T1 += T
        T2 += pow(T, 2)
        T3 += pow(T, 3)
        T4 += pow(T, 4)
        T5 += pow(T, 5)
        T6 += pow(T, 6)

        RT1 += R * T / 1.0
        RT2 += R * T * T / 2.0
        RT3 += R * T * T * T / 3.0

    A = np.array([
        [T6 / 9.0, T5 / 6.0, T4 / 3.0],
        [T5 / 6.0, T4 / 4.0, T3 / 2.0],
        [T4 / 3.0, T3 / 2.0, T2 / 1.0]
    ])

    #print("A: ", A)

    b = np.array([RT3, RT2, RT1])

    x = np.linalg.solve(A, b)
    if (normalize):
        x = x / x[2]

    # test inversion
    # print("b: ", b)
    # print("Ax: ", np.matmul(A, x) )

    return x

def calculateRotationSpeedFromCoeff(c, Tary):
    w = []
    for T in Tary:
        w.append(pow(T, 2) * c[0] + T * c[1] + c[2])
    return w

def calculateTotalRotationFromCoeff(c, Tary):
    R = []
    for T in Tary:
        R.append((float)(pow(T,3) * c[0] / 3.0 + pow(T,2) * c[1] / 2.0 + T *c[2]))
    return R

def calculateRegressionRmsError(c, RTlst):
    RList = [row[0] for row in RTlst]
    TList = [row[1] for row in RTlst]
    Rc = calculateTotalRotationFromCoeff(c, TList)
    N = len(TList)
    sum = 0
    for i in range(N):
        e = RList[i] - Rc[i]
        sum = e * e
    rms = np.sqrt(sum) / N
    return rms;
         
def createTestDataSet(a, b, w_0, tMax, N):
    RT = []
    deltaT = tMax / N
    for i in range(N):
        T = i * deltaT
        R = (1.0/3.0 * a * pow(T,3) + 0.5 * b * pow(T,2) + T) * w_0
        RT.append([R, T])
    return RT

if __name__ == "__main__":
    # The return value of main() is passed to sys.exit()
    # to provide a proper exit code to the operating system.
    sys.exit(main())