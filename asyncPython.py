from multiprocessing import *
import os
import numpy as np


SI=['k','p','q','s','u','v','x','g','n','r','aa','ab']
SD=['a','b','f','h','l','ma','c','d','e','i','j','o','t','w','y','z']

EI=['ac','ad']
ED=['ae']


def normalize_params(sequence_list):
    print("starting normalisation of read sequences")
    pool = Pool()
    param_list = pool.map(calculateParameters, sequence_list)
    pool.close()
    return param_list

def energyStructParamsMP(normalised_params_list):

    nml = len(normalised_params_list)
    print("Starting combining Energy and Struct")
    pool = Pool()
    SIParams_all_seq = pool.starmap(combineStructEnergyParams,[(SI,list(normalised_params_list[seq].items())) for seq in range(nml)])
    pool.close()
    pool.join()

    pool = Pool()
    SDParams_all_seq = pool.starmap(combineStructEnergyParams,[(SD,list(normalised_params_list[seq].items())) for seq in range(nml)])
    pool.close()
    pool.join()

    pool = Pool()
    EIparams_all_seq = pool.starmap(combineStructEnergyParams,[(EI,list(normalised_params_list[seq].items())) for seq in range(nml)])
    pool.close()
    pool.join()

    pool = Pool()
    EDParams_all_seq = pool.starmap(combineStructEnergyParams,[(ED,list(normalised_params_list[seq].items())) for seq in range(nml)])
    pool.close()
    pool.join()

    combined_params_map = dict(zip(['SIParams_all_seq','SDParams_all_seq', 'EIparams_all_seq', 'EDParams_all_seq'],[SIParams_all_seq,SDParams_all_seq, EIparams_all_seq, EDParams_all_seq]))
    return transformStructEnerMap(combined_params_map)

def assign_params(param_map,a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab, ac, ad, ae):
    param_map['a'].append(a)
    param_map['b'].append(b)
    param_map['c'].append(c)
    param_map['d'].append(d)
    param_map['e'].append(e)
    param_map['f'].append(f)
    param_map['g'].append(g)
    param_map['h'].append(h)
    param_map['i'].append(i)
    param_map['j'].append(j)
    param_map['k'].append(k)
    param_map['l'].append(l)
    param_map['ma'].append(ma)
    param_map['n'].append(n)
    param_map['o'].append(o)
    param_map['p'].append(p)
    param_map['q'].append(q)
    param_map['r'].append(r)
    param_map['s'].append(s)
    param_map['t'].append(t)
    param_map['u'].append(u)
    param_map['v'].append(v)
    param_map['w'].append(w)
    param_map['x'].append(x)
    param_map['y'].append(y)
    param_map['z'].append(z)
    param_map['aa'].append(aa)
    param_map['ab'].append(ab)
    param_map['ac'].append(ac)
    param_map['ad'].append(ad)
    param_map['ae'].append(ae)
    return param_map

def calculateMovingAverages(param_map):
    moving_win_size = 25
    moving_param_map = {}
    for k,v in param_map.items():
        arr = v
        moving_param_map[k] = []
        for i in range(0, len(arr)-moving_win_size+1):
            this_window = arr[i : i+moving_win_size]
            window_avg = sum(this_window)/moving_win_size
            moving_param_map[k].append(window_avg)
    return normalizeMovingAverages(moving_param_map)

def normalizeMovingAverages(moving_param_map):
    normalized_map = {}
    for k in moving_param_map.keys():
        arr = moving_param_map[k]
        minArr = min(arr)
        range = max(arr)-minArr
        normalized_map[k] = []
        for i in arr:
            norm_val = (i-minArr)/range
            normalized_map[k].append(norm_val)
    return normalized_map

def calculateParameters(b_arr2):
    param_map = {'a':[],'b':[],'c':[],'d':[],'e':[],'f':[],'g':[],'h':[],'i':[],'j':[],'k':[],'l':[],'ma':[],'n':[],'o':[],'p':[],'q':[],'r':[],'s':[],'t':[],'u':[],'v':[],'w':[],'x':[],'y':[],'z':[],'aa':[],'ab':[],'ac':[],'ad':[],'ae':[]}
    noofbases = len(b_arr2)
    a=b=c=d=e=f=g=h=i=j=k=l=ma=n=o=p=q=r=s=t=u=v=w=x=y=z=aa=ab=ac=ad=ae = 0
    if noofbases == 0:
        return
    for m in range(noofbases-1):
        if b_arr2[m] == 'A' and b_arr2[m+1] == 'A' or b_arr2[m] == 'T' and b_arr2[m+1] == 'T':
            a = -0.24586
            b = 0.004276
            c = -0.61382
            d = 0.596711
            e = 1.421711
            f = 0.128291
            g = -0.17956
            h = 0.120633
            i = -1.48165
            j = -15.3082
            k = 3.158861
            l = -0.11
            ma = -0.2
            n = 3.25
            o = 0.63
            p = -0.08
            q = 35.67
            r = 3.27
            s = 35.56
            t = -55.5495
            u = 51.04402
            v = 50.77228
            w = 128.4076
            x = -5.92446
            y = -96.656
            z = -110.254
            aa = 118.7842
            ab = 38.99565
            ac = -5.44
            ad = -26.71
            ae = -171.84
            assign_params(param_map,a,b,c,d,e,f,g,h,i,j,k,l,ma,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae)
        elif (b_arr2[m] ==  'G' and b_arr2[m+1] == 'G' or (b_arr2[m] ==  'C' and b_arr2[m+1] == 'C')):
            a = 0.584605
            b = -0.07763
            c = 9.909211
            d = 3.127632
            e = 2.585526
            f = 1.2239
            g = -1.1225
            h = 1.6701
            i = 8.842
            j = 5.707
            k = 5.719
            l = 0.38
            ma = 0.71
            n = 2.96
            o = -1.62
            p = -2.55
            q = 29.49
            r = 2.73
            s = 31.57
            t = -38.322
            u = 39.807
            v = 44.997
            w = 129.998
            x = -56.716
            y = -61.65
            z = -103.913
            aa = 96.172
            ab = 37.361
            ac = -8.48
            ad = -26.28
            ae = -166.76
            assign_params(param_map,a,b,c,d,e,f,g,h,i,j,k,l,ma,n,o,p,q,r,s,t,u,v, w, x, y,z,aa,ab,ac,ad,ae)
        elif ((b_arr2[m].upper() == 'T') and (b_arr2[m+1].upper() == 'A')):
            a = 0.101842
            b = -0.28579
            c = 0.547368
            d = -0.76579
            e = 1.471053
            f = 0.026053
            g = -0.17737
            h = 0.084737
            i = -0.32632
            j = -14.0368
            k = 3.305263
            l = -0.05
            ma = 0.37
            n = 3.39
            o = -2.71
            p = 1.74
            q = 32.05
            r = 3.41
            s = 32
            t = -50.82
            u = 29.23
            v = 53.74
            w = 127.846
            x = -37.772
            y = -61.97
            z = -106.74
            aa = 119.92
            ab = 37.824
            ac = -5.83
            ad = -26.9
            ae = -174.35
            assign_params(param_map,a,b,c,d,e,f,g,h,i,j,k,l,ma,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae)
        elif(b_arr2[m].upper()=='C' and b_arr2[m+1].upper()=='G'):
            a = 0.70475
            b = 0.060125
            c = 2.910625
            d = 0.303125
            e = 1.48875
            f = 0.355364
            g = -0.35755
            h = 0.561455
            i = 2.413182
            j = -3.27955
            k = -0.84909
            l = 0.24
            ma = 0.68
            n = 3.33
            o = 1.75
            p = 4.29
            q = 37.38
            r = 3.33
            s = 39.1
            t = -43.178
            u = 57.555
            v = 51.502
            w = 132.616
            x = -67.013
            y = -62.643
            z = -106.586
            aa = 103.3
            ab = 37.951
            ac = -8.05
            ad = -27.93
            ae = -176.88
            assign_params(param_map,a,b,c,d,e,f,g,h,i,j,k,l,ma,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae)
        elif b_arr2[m].upper() == 'G' and b_arr2[m+1].upper() == 'C':
            a = 0.593382
            b = 0.082647
            c = 0.364706
            d = 2.507353
            e = 1.244118
            f = 0.021286
            g = -0.13686
            h = 0.284
            i = -0.98286
            j = -10.4429
            k = -1.11571
            l = -0.28
            ma = 0.17
            n = 3.35
            o = -2.68
            p = -5.28
            q = 31.09
            r = 3.29
            s = 30.56
            t = -48.267
            u = 36.968
            v = 36.646
            w = 129.951
            x = -70.502
            y = -63.978
            z = -109.979
            aa = 84.877
            ab = 38.594
            ac = -8.72
            ad = -28.13
            ae = -165.58
            assign_params(param_map,a,b,c,d,e,f,g,h,i,j,k,l,ma,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae)
        elif (b_arr2[m].upper() == 'A' and b_arr2[m+1].upper() == 'T'):
            a = -0.17819
            b = -0.02597
            c = -0.25972
            d = 0.461111
            e = 1.2125
            f = 0.01375
            g = -0.13375
            h = 0.147222
            i = 0.279167
            j = -16.5361
            k = 4.109722
            l = -0.06
            ma = -0.44
            n = 3.31
            o = -0.9
            p = -2.62
            q = 33.42
            r = 3.33
            s = 33.68
            t = -51.913
            u = 59.981
            v = 47.175
            w = 125.847
            x = -31.23
            y = -93.609
            z = -113.738
            aa = 107.8554
            ab = 40.016
            ac = -5.35
            ad = -27.2
            ae = -173.7
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif(b_arr2[m]=='A' and b_arr2[m+1] == 'C') or (b_arr2[m] == 'G' and b_arr2[m+1] == 'T'):
            a = -0.03125
            b = -0.0505
            c = 0.68
            d = 2.51
            e = 1.25
            f = 0.055
            g = -0.12481
            h = 0.572885
            i = -5.57885
            j = -9.83846
            k = 2.109615
            l = -0.12
            ma = 0.73
            n = 3.66
            o = -4.75
            p = -0.6
            q = 34.79
            r = 2.97
            s = 37.6
            t = -45.4742
            u = 42.118
            v = 51.8606
            w = 126.878
            x = -34.3212
            y = -55.7666
            z = -116.5469
            aa = 106.334
            ab = 38.65
            ac = -7.14
            ad = -27.73
            ae = -171.11
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif (b_arr2[m] ==  'C' and b_arr2[m+1] == 'A') or (b_arr2[m] ==  'T' and b_arr2[m+1] == 'G'):
            a = 0.311563
            b = -0.0675
            c = 3.004688
            d = 0.175
            e = 1.054688
            f = -0.00182
            g = -0.14879
            h = 0.156212
            i = 0.075758
            j = -10.8242
            k = 2.072727
            l = -0.07
            ma = 0.52
            n = 3.29
            o = -3.27
            p = 1.54
            q = 36.74
            r = 3.21
            s = 37.25
            t = -52.023
            u = 41.3456
            v = 35.389
            w = 131.452
            x = -75.781
            y = -60.38
            z = -91.5032
            aa = 73.933
            ab = 38.744
            ac = -7.01
            ad = -27.15
            ae = -179.01
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w,x, y, z, aa, ab,
                          ac, ad, ae)
        elif (b_arr2[m] ==  'A' and b_arr2[m+1] == 'G') or (b_arr2[m] ==  'C' and b_arr2[m+1] == 'T'):
            a = 0.485909
            b = -0.02864
            c = 1.825
            d = 3.722727
            e = 1.652273
            f = -0.06364
            g = -0.16818
            h = -0.14614
            i = -0.85
            j = -11.7341
            k = -3.66818
            l = -0.27
            ma = 0.21
            n = 3.02
            o = 3.69
            p = -7.27
            q = 31
            r = 3.06
            s = 19.06
            t = -36.214
            u = 36.155
            v = 43.87
            w = 129.955
            x = -43.587
            y = -34.57
            z = -88.341
            aa = 73.08
            ab = 38.941
            ac = -6.27
            ad = -26.89
            ae = -174.93
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif (b_arr2[m] ==  'G' and b_arr2[m+1] == 'A') or (b_arr2[m] ==  'T' and b_arr2[m+1] == 'C'):
            a = -0.15833
            b = 0.08
            c = 0.525926
            d = 0.751852
            e = 1.177778
            f = 0.049815
            g = -0.19926
            h = 0.21463
            i = -1.16852
            j = -13.1352
            k = 0.942593
            l = 0.04
            ma = -0.2
            n = 3.32
            o = 2.21
            p = 2.69
            q = 37.74
            r = 3.33
            s = 37.81
            t = -42.611
            u = 54.835
            v = 53.374
            w = 131.432
            x = -8.081
            y = -94.004
            z = -110.983
            aa = 118.188
            ab = 37.945
            ac = -7.8
            ad = -26.78
            ae = -167.6
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif (b_arr2[m]=='A' and b_arr2[m] == 'N'):
            a = 0.217
            b = -0.031
            c = 1.889
            d = 1.339
            e = 1.456
            f = 0.181
            g = -0.275
            h = 0.367
            i = 0.122
            j = -9.943
            k = 1.578
            l = -0.030
            ma = 0.255
            n = 3.288
            o = -0.765
            p = -0.814
            q = 33.937
            r = 3.193
            s = 33.419
            t = -46.437
            u = 44.904
            v = 46.933
            w = 129.438
            x = -43.093
            y = -68.523
            z = -105.858
            aa = 100.244
            ab = 38.502
            ac = -7.009
            ad = -27.170
            ae = -172.176
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif b_arr2[m] == 'T' and b_arr2[m+1] == "N":
            a = 0.217
            b = -0.031
            c = 1.889
            d = 1.339
            e = 1.456
            f = 0.181
            g = -0.275
            h = 0.367
            i = 0.122
            j = -9.943
            k = 1.578
            l = -0.030
            ma = 0.255
            n = 3.288
            o = -0.765
            p = -0.814
            q = 33.937
            r = 3.193
            s = 33.419
            t = -46.437
            u = 44.904
            v = 46.933
            w = 129.438
            x = -43.093
            y = -68.523
            z = -105.858
            aa = 100.244
            ab = 38.502
            ac = -7.009
            ad = -27.170
            ae = -172.176
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif (b_arr2[m]== 'G' and b_arr2[m+1] == 'N'):
            a = 0.217
            b = -0.031
            c = 1.889
            d = 1.339
            e = 1.456
            f = 0.181
            g = -0.275
            h = 0.367
            i = 0.122
            j = -9.943
            k = 1.578
            l = -0.030
            ma = 0.255
            n = 3.288
            o = -0.765
            p = -0.814
            q = 33.937
            r = 3.193
            s = 33.419
            t = -46.437
            u = 44.904
            v = 46.933
            w = 129.438
            x = -43.093
            y = -68.523
            z = -105.858
            aa = 100.244
            ab = 38.502
            ac = -7.009
            ad = -27.170
            ae = -172.176
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif(b_arr2[m]=='C' and b_arr2[m+1]=='N'):
            a = 0.217
            b = -0.031
            c = 1.889
            d = 1.339
            e = 1.456
            f = 0.181
            g = -0.275
            h = 0.367
            i = 0.122
            j = -9.943
            k = 1.578
            l = -0.030
            ma = 0.255
            n = 3.288
            o = -0.765
            p = -0.814
            q = 33.937
            r = 3.193
            s = 33.419
            t = -46.437
            u = 44.904
            v = 46.933
            w = 129.438
            x = -43.093
            y = -68.523
            z = -105.858
            aa = 100.244
            ab = 38.502
            ac = -7.009
            ad = -27.170
            ae = -172.176

            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif(b_arr2[m] == 'N' and b_arr2[m+1] == 'A'):
            a = 0.217
            b = -0.031
            c = 1.889
            d = 1.339
            e = 1.456
            f = 0.181
            g = -0.275
            h = 0.367
            i = 0.122
            j = -9.943
            k = 1.578
            l = -0.030
            ma = 0.255
            n = 3.288
            o = -0.765
            p = -0.814
            q = 33.937
            r = 3.193
            s = 33.419
            t = -46.437
            u = 44.904
            v = 46.933
            w = 129.438
            x = -43.093
            y = -68.523
            z = -105.858
            aa = 100.244
            ab = 38.502
            ac = -7.009
            ad = -27.170
            ae = -172.176
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif (b_arr2[m] ==  'N' and b_arr2[m+1] == 'T'):
            a = 0.217
            b = -0.031
            c = 1.889
            d = 1.339
            e = 1.456
            f = 0.181
            g = -0.275
            h = 0.367
            i = 0.122
            j = -9.943
            k = 1.578
            l = -0.030
            ma = 0.255
            n = 3.288
            o = -0.765
            p = -0.814
            q = 33.937
            r = 3.193
            s = 33.419
            t = -46.437
            u = 44.904
            v = 46.933
            w = 129.438
            x = -43.093
            y = -68.523
            z = -105.858
            aa = 100.244
            ab = 38.502
            ac = -7.009
            ad = -27.170
            ae = -172.176
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif (b_arr2[m] ==  'N' and b_arr2[m+1] == 'G'):
            a = 0.217
            b = -0.031
            c = 1.889
            d = 1.339
            e = 1.456
            f = 0.181
            g = -0.275
            h = 0.367
            i = 0.122
            j = -9.943
            k = 1.578
            l = -0.030
            ma = 0.255
            n = 3.288
            o = -0.765
            p = -0.814
            q = 33.937
            r = 3.193
            s = 33.419
            t = -46.437
            u = 44.904
            v = 46.933
            w = 129.438
            x = -43.093
            y = -68.523
            z = -105.858
            aa = 100.244
            ab = 38.502
            ac = -7.009
            ad = -27.170
            ae = -172.176
            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                          ac, ad, ae)
        elif (b_arr2[m] ==  'N' and b_arr2[m+1] == 'C'):
            a = 0.217
            b = -0.031
            c = 1.889
            d = 1.339
            e = 1.456
            f = 0.181
            g = -0.275
            h = 0.367
            i = 0.122
            j = -9.943
            k = 1.578
            l = -0.030
            ma = 0.255
            n = 3.288
            o = -0.765
            p = -0.814
            q = 33.937
            r = 3.193
            s = 33.419
            t = -46.437
            u = 44.904
            v = 46.933
            w = 129.438
            x = -43.093
            y = -68.523
            z = -105.858
            aa = 100.244
            ab = 38.502
            ac = -7.009
            ad = -27.170
            ae = -172.176

            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,
                      ac, ad, ae)
        elif (b_arr2[m]=='N' and b_arr2[m+1] == 'N'):
            a = 0.217
            b = -0.031
            c = 1.889
            d = 1.339
            e = 1.456
            f = 0.181
            g = -0.275
            h = 0.367
            i = 0.122
            j = -9.943
            k = 1.578
            l = -0.030
            ma = 0.255
            n = 3.288
            o = -0.765
            p = -0.814
            q = 33.937
            r = 3.193
            s = 33.419
            t = -46.437
            u = 44.904
            v = 46.933
            w = 129.438
            x = -43.093
            y = -68.523
            z = -105.858
            aa = 100.244
            ab = 38.502
            ac = -7.009
            ad = -27.170
            ae = -172.176

            assign_params(param_map, a, b, c, d, e, f, g, h, i, j, k, l, ma, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab,ac, ad, ae)
    return calculateMovingAverages(param_map)

def combineStructEnergyParams(array,normalized_list_tuples):
    normalized_map = dict(normalized_list_tuples)
    map = np.zeros(976)
    for k in array:
        arr = normalized_map[k]
        for i in range(len(arr)):
            map[i]+=arr[i]
    return map

def transformStructEnerMap(struct_ener_map):
    transformed_map = {}
    for k in struct_ener_map.keys():
        values_of_all_seq_per_param = struct_ener_map[k]
        for i in range(len(values_of_all_seq_per_param)):
            try:
                transformed_map[i]
            except KeyError:
                transformed_map[i] = {}
            transformed_map[i][k] = values_of_all_seq_per_param[i]
    return transformed_map