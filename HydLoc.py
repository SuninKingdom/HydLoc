from tkinter import *
from tkinter import scrolledtext
import tkinter.font as tf
import pickle
import numpy as np

# *******************************************************************************************************************
# some constants used in the script
# 1:P, 2:K
F1 = np.array([0.067986, 0.0135, 0.049038, 0.059339, 0.029195, 0.157604, 0.015533,
               0.03589, 0.049065, 0.072811, 0.015777, 0.030062, 0.113039, 0.040499,
               0.049607, 0.064652, 0.050095, 0.059013, 0.007184, 0.020114])

F2 = np.array([0.069633, 0.010413, 0.049122, 0.05943, 0.027138, 0.190596, 0.011991,
               0.032713, 0.048122, 0.063322, 0.013148, 0.026822, 0.128326, 0.038393,
               0.054328, 0.054644, 0.043442, 0.055906, 0.005049, 0.017461])

PC = np.array([[0.62,0.29,-0.9,-0.74,1.19,0.48,-0.4,1.38,-1.5,1.06,0.64,-0.78,0.12,-0.85,-2.53,-0.18,-0.05,1.08,0.81,0.26],
      [-0.5,-1,3,3,-2.5,0,-0.5,-1.8,3,-1.8,-1.3,2,0,0.2,3,0.3,-0.4,-1.5,-3.4,-2.3],
      [27.5,44.6,40,62,115.5,0,79,93.5,100,93.5,94.1,58.7,41.9,80.7,105,29.3,51.3,71.5,145.5,117.3],
      [8.1,5.5,13,12.3,5.2,9,10.4,5.2,11.3,4.9,5.7,11.6,8,10.5,10.5,9.2,8.6,5.9,5.4,6.2],
      [0.046,0.128,0.105,0.151,0.29,0,0.23,0.186,0.219,0.186,0.221,0.134,0.131,0.18,0.291,0.062,0.108,0.14,0.409,0.298],
      [1.181,1.461,1.587,1.862,2.228,0.881,2.025,1.81,2.258,1.931,2.034,1.655,1.468,1.932,2.56,1.298,1.525,1.645,2.663,2.368],
      [0.007187,-0.03661,-0.02382,0.006802,0.037552,0.179052,-0.01069,0.021631,0.017708,0.051672,0.002683,0.005392,0.239531,
       0.049211,0.043587,0.004627,0.003352,0.057004,0.037977,0.023599]])

P = np.array([['RKEDQN','GASTPHY', 'CVLIMFW'],['GASCTPD', 'NVEQIL', 'MHKFRYW'],
              ['LIFWCMVY', 'PATGS', 'HQRKNED'], ['GASDT', 'CPNVEQIL', 'KMHFRYW']])

aa = 'ACDEFGHIKLMNPQRSTVWY'

with open('clf1', 'rb') as m1:
    clf1 = pickle.load(m1)

with open('clf2', 'rb') as m2:
    clf2 = pickle.load(m2)

RT = 0

AC = []
length = []            # length stores the length of protein
result1 = []           # result stores site, label and probability information
result2 = []
frag1 = []             # frag stores the information of peptide samples
frag2 = []

# *******************************************************************************************************************
# Here are some functions:
def show_menu(e):
    global m
    w = e.widget
    m.entryconfigure("Copy",
    command=lambda: w.event_generate("<<Copy>>"))
    m.entryconfigure("Cut",
    command=lambda: w.event_generate("<<Cut>>"))
    m.entryconfigure("Paste",
    command=lambda: w.event_generate("<<Paste>>"))
    m.tk.call("tk_popup", m, e.x_root, e.y_root)

def find_ss(s):
    l = len(s)
    site_P = []
    site_K = []
    pep_P = []
    pep_K = []
    for k in range(len(s)):
        if 8<k<l-9 and s[k] == 'P':
            site_P.append(k)
            pep_P.append(s[k-9: k+10])
        elif 13<k<l-14 and s[k] == 'K':
            site_K.append(k)
            pep_K.append(s[k-14: k+15])
    return site_P, pep_P, site_K, pep_K

def repla_ambi(s):
    global aa
    for i2 in range(len(s)):
        if s[i2] not in aa:
            ind1 = (i2+21)%20
            s[i2] = aa[ind1]
    return s

def cal_AACre(s, t=1):
    global F1, F2
    aac = np.zeros((1, 20))
    for i in range(20):
        n = s.count(aa[i])
        f = n/len(s)
        aac[0][i] = f

    if t == 1:
        F = F1
    else:
        F = F2
    r1 = aac/F
    r2 = np.ones((20, ))
    for j in range(20):
        if r1[0][j] != 0:
            r2[j] = r1[0][j]

    lg = np.log2(r2)
    p = aac * lg
    re = p.sum()
    return aac, re

def cal_PSP(s):
    global aa, PC
    psp1 = np.zeros((7, len(s)))
    for i in range(len(s)):
        ind = aa.index(s[i])
        psp1[:, i] = PC[:, ind]
    psp2 = psp1.reshape((1, 7*len(s)))
    return psp2

def cal_CTD(s):
    '''
    :param s: the amino acid sequence, the length of the sequence must bigger than 1
    :return: composition(c), transition(t), distribution(d)
    '''
    global P # P were four properties used in ctd, each property has 3 types
    l = len(s)
    c1 = np.zeros((4, 3))
    t1 = np.zeros((4, 3))
    d1 = np.zeros((4, 15))
    s2 = np.zeros((4, len(s)), dtype=str) # s2 was a sequence of symbols representing different properties
    # first, change the a.a. sequence to symbols
    for i in range(l):
        r = s[i]
        for j in range(4):
            p = P[j, :]
            for k in range(3):
                if r in p[k]:
                    s2[j, i] = k
    # then, find the indexes of a specific type in a specific property
    for ii in range(4):
        s3 = ''.join(s2[ii, :])
        ind1 = []  # the index of type I of a specific property
        ind2 = []
        ind3 = []
        ind = []   # the collection of ind1, ind2, ind3
        for kk in range(l):
            if s3[kk] == '0':
                ind1.append(kk)
            elif s3[kk] == '1':
                ind2.append(kk)
            else:
                ind3.append(kk)
        ind.insert(0, ind1)
        ind.insert(1, ind2)
        ind.insert(2, ind3)
    # finally, calculate c,t,d
        for jj in range(3):
            ind4 = ind[jj]
            n = len(ind4)
            c1[ii, jj] = n/l
            if n != 0:
                d1[ii, 5*jj] = (ind4[0]+1)/l
                m1 = round(n*0.25 + 1e-8) # because want round(0.5) = 1
                if m1 != 0:
                    d1[ii, 1+5*jj] = (ind4[m1-1]+1)/l
                m2 = round(n*0.5 + 1e-8)-1
                d1[ii, 2+5*jj] = (ind4[m2]+1)/l
                d1[ii, 3+5*jj] = (ind4[round(n*0.75+1e-8)-1]+1)/l
                d1[ii, 4+5*jj] = (ind4[-1] + 1) / l
            del ind4, n
        t1[ii, 0] = (s3.count('01') + s3.count('10')) / (l-1)
        t1[ii, 1] = (s3.count('02')+s3.count('20')) / (l-1)
        t1[ii, 2] = (s3.count('12')+s3.count('21')) / (l-1)
        del s3
    c = c1.reshape((1, 12))
    t = t1.reshape((1, 12))
    d = d1.reshape((1, 60))
    ctd = np.c_[c, t, d]
    return c, t, d, ctd

def ChooseR():
    '''
    :return: residue type, as a parameter of pred()
    '''
    global RT
    RT = v.get()
    print(RT)
    return RT

def pred(Rt=0):
    global RT, tp, AC, length, frag1, frag2, result1, result2
    Rt = RT
    c1 = t1.get(0.0, "end")
    print(len(c1))
    c2 = c1.split('>')
    print(c2)
    c2.pop(0)
    print(c2)
    AC = []
    length = []
    frag1 = []
    frag2= []
    result1 = []
    result2 = []
    seq = []
    for i in range(len(c2)):
        c3 = c2[i]
        c4 = c3.split('\n')
        AC.append(c4[0])
        c4.pop(0)
        c5 = ''.join(c4)
        seq.append(c5)

    length = [len(le) for le in seq]

    # make prediction and plot the results according to the selected residue
    # set some font type first

    wlist = root.winfo_children()
    wclass = [wid.winfo_class() for wid in wlist]

    if 'Toplevel' in wclass:
        tp.destroy()

    tp = Toplevel()
    tp.configure(bg='white')

    # make prediction and save the result first
    for j in range(len(AC)):
        [si1, pep1, si2, pep2] = find_ss(seq[j])
        n1 = len(si1)
        n2 = len(si2)
        pep11 = [repla_ambi(p1) for p1 in pep1]
        pep22 = [repla_ambi(p2) for p2 in pep2]
        res1 = np.zeros((n1, 3))
        res2 = np.zeros((n2, 3))
        if Rt == 0:              # make prediction of P residue
            if n1 != 0:
                for jj in range(n1):
                    pep111 = pep11[jj]
                    res1[jj, 0] = si1[jj]
                    [aac1, re1] = cal_AACre(pep111)
                    psp1 = cal_PSP(pep111)
                    [_, _, _, ctd1] = cal_CTD(pep111)
                    fea1 = np.c_[aac1, re1, psp1, ctd1]
                    res1[jj, 1] = clf1.predict(fea1)
                    print(fea1)
                    res1[jj, 2] = clf1.predict_proba(fea1)[0, 1]
        result1.append(res1)
        frag1.append(pep11)

        if Rt == 1:              # make prediction of K residue
            if n2 != 0:
                for jj in range(n2):
                    pep222 = pep22[jj]
                    res2[jj, 0] = si2[jj]
                    [aac2, re2] = cal_AACre(pep222, 2)
                    psp2 = cal_PSP(pep222)
                    [_, _, _, ctd2] = cal_CTD(pep222)
                    fea2 = np.c_[aac2, re2, psp2, ctd2]
                    res2[jj, 1] = clf2.predict(fea2)
                    print(fea2)
                    res2[jj, 2] = clf2.predict_proba(fea2)[0, 1]
        result2.append(res2)
        frag2.append(pep22)

        if Rt == 2:              # make prediction of P&K residue
            if n1 != 0:
                for jj in range(n1):
                    pep111 = pep11[jj]
                    res1[jj, 0] = si1[jj]
                    [aac1, re1] = cal_AACre(pep111)
                    psp1 = cal_PSP(pep111)
                    [_, _, _, ctd1] = cal_CTD(pep111)
                    fea1 = np.c_[aac1, re1, psp1, ctd1]
                    res1[jj, 1] = clf1.predict(fea1)
                    print(fea1)
                    res1[jj, 2] = clf1.predict_proba(fea1)[0, 1]

            if n2 != 0:
                for jj in range(n2):
                    pep222 = pep22[jj]
                    res2[jj, 0] = si2[jj]
                    [aac2, re2] = cal_AACre(pep222)
                    psp2 = cal_PSP(pep222)
                    [_, _, _, ctd2] = cal_CTD(pep222,2)
                    fea2 = np.c_[aac2, re2, psp2, ctd2]
                    res2[jj, 1] = clf2.predict(fea2)
                    print(fea2)
                    res2[jj, 2] = clf2.predict_proba(fea2)[0, 1]

        result1.append(res1)
        frag1.append(pep11)
        result2.append(res2)
        frag2.append(pep22)

        # project the result on the toplevel window
        Label(tp, text='Prediction  result', font=('Verdana', 17, 'bold')).grid(row=0, column=0)
        if Rt == 0:
            r1 = 0
            r2 = 0
            r3 = 0
            r4 = 0
            for k in range(len(AC)):
                resu1 = result1[k]
                n11 = resu1.shape[0]
                pep_1 = frag1[k]
                if n11 != 0:
                    ac = AC[k]
                    l = length[k]
                    e00 = Entry(tp, borderwidth=0, foreground='Blue',
                                font=('Arial', 16, 'bold'))
                    e01 = Entry(tp, borderwidth=0, font=('Arial', 12))
                    e1 = Entry(tp, borderwidth=0, font=('Vani', 12, 'bold'))
                    e20 = Entry(tp, borderwidth=0, font=('Cooper Black', 11))
                    e21 = Entry(tp, borderwidth=0, font=('Cooper Black', 11))
                    e22 = Entry(tp, borderwidth=0, font=('Cooper Black', 11))
                    e23 = Entry(tp, borderwidth=0, font=('Cooper Black', 11))
                    e00.grid(row=r1, column=0)
                    e01.grid(row=r1, column=1)
                    r1 = r1 + 1
                    e1.grid(row=r1, column=0)
                    r1 = r1 + 1
                    e20.grid(row=r1, column=0)
                    e21.grid(row=r1, column=1)
                    e22.grid(row=r1, column=2)
                    e23.grid(row=r1, column=3)
                    r1 = r1 + 1
                    e00.insert(0, '$' + ac)
                    e01.insert(0, '(The length is ' + str(l) + ' AAs)')
                    e1.insert(0, 'Hydroxyproline')
                    e20.insert(0, 'site')
                    e21.insert(0, 'Sequence')
                    e22.insert(0, 'Probability')
                    e23.insert(0, 'Yes/No')
                    for kk in range(n11):
                        site1 = resu1[kk, 0]
                        label1 = resu1[kk, 1]
                        proba1 = resu1[kk, 2]
                        pep_11 = pep_1[kk]
                        if label1 == 1:
                            e_0 = Entry(tp, borderwidth=0, foreground='Red',
                                        font=('Times', 12))
                            e_1 = Entry(tp, borderwidth=0, font=('Times', 12))
                            e_2 = Entry(tp, borderwidth=0, font=('Times', 12))
                            e_3 = Entry(tp, borderwidth=0, foreground='Red',
                                        font=('Times', 12))
                            e_0.grid(row=r1, column=0)
                            e_1.grid(row=r1, column=1)
                            e_2.grid(row=r1, column=2)
                            e_3.grid(row=r1, column=3)
                            e_0.insert(0, str(int(site1+1)))
                            e_1.insert(0, pep_11)
                            e_2.insert(0, str(round(proba1, 4)))
                            e_3.insert(0, 'Yes')
                            r1 = r1 + 1
                        else:
                            e_0 = Entry(tp, borderwidth=0, font=('Times', 12))
                            e_1 = Entry(tp, borderwidth=0, font=('Times', 12))
                            e_2 = Entry(tp, borderwidth=0, font=('Times', 12))
                            e_3 = Entry(tp, borderwidth=0, font=('Times', 12))
                            e_0.grid(row=r1, column=0)
                            e_1.grid(row=r1, column=1)
                            e_2.grid(row=r1, column=2)
                            e_3.grid(row=r1, column=3)
                            e_0.insert(0, str(int(site1+1)))
                            e_1.insert(0, pep_11)
                            e_2.insert(0, str(round(proba1, 4)))
                            e_3.insert(0, 'No')
                            r1 = r1 + 1

    return AC, length, frag1, frag2, result1, result2

# *******************************************************************************************************************
# Here the interface design

root = Tk()
root.configure(bg='black')
m = Menu(root, tearoff=0)
m.add_command(label="Copy")
m.add_command(label="Cut")
m.add_command(label="Paste")

t1 = scrolledtext.ScrolledText(root, width=80, height=18, bg='#F8F8FF')
t1.grid(row=0, column=0, sticky=N)

v = IntVar()
rtv = [("P", 0), ("K", 1), ("P&K", 2)]  # residue types and corresponding values

for rt, rv in rtv:
    rb = Radiobutton(root, text=rt, value=rv, command=ChooseR, variable=v, font=('Vani', 10, 'bold'))
    rb.grid(row=rv+2, column=0, sticky=W)

bt = Button(root, text='Submit', width=10, height=2, command=pred,
            bg='#C0C0C0', font=('Verdana', 10, 'bold'))
bt.grid(row=5, column=0)

lb1 = Label(root, text='Choose residues: ',
            font=('David', 13, 'bold')).grid(row=1, column=0, sticky=W)

t1.bind_class("Text", "<Button-3><ButtonRelease-3>", show_menu)

root.mainloop()
