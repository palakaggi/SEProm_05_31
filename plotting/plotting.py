import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path


def plotting(parameter_map_tss, parameter_map_cds, name):
    len_seq_param = len(parameter_map_tss['normalized_params_map'][0]['a'])

    a = [0 for i in range(len_seq_param)]
    b = [0 for i in range(len_seq_param)]
    c = [0 for i in range(len_seq_param)]
    d = [0 for i in range(len_seq_param)]
    e = [0 for i in range(len_seq_param)]
    f = [0 for i in range(len_seq_param)]
    g = [0 for i in range(len_seq_param)]
    h = [0 for i in range(len_seq_param)]
    i = [0 for m in range(len_seq_param)]
    j = [0 for i in range(len_seq_param)]
    k = [0 for i in range(len_seq_param)]
    l = [0 for i in range(len_seq_param)]
    m = [0 for i in range(len_seq_param)]
    n = [0 for i in range(len_seq_param)]
    o = [0 for i in range(len_seq_param)]
    p = [0 for i in range(len_seq_param)]
    q = [0 for i in range(len_seq_param)]
    r = [0 for i in range(len_seq_param)]
    s = [0 for i in range(len_seq_param)]
    t = [0 for i in range(len_seq_param)]
    u = [0 for i in range(len_seq_param)]
    v = [0 for i in range(len_seq_param)]
    w = [0 for i in range(len_seq_param)]
    x = [0 for i in range(len_seq_param)]
    y = [0 for i in range(len_seq_param)]
    z = [0 for i in range(len_seq_param)]
    aa = [0 for i in range(len_seq_param)]
    ab = [0 for i in range(len_seq_param)]
    ac = [0 for i in range(len_seq_param)]
    ad = [0 for i in range(len_seq_param)]
    ae = [0 for i in range(len_seq_param)]

    a_cds = [0 for i in range(len_seq_param)]
    b_cds = [0 for i in range(len_seq_param)]
    c_cds = [0 for i in range(len_seq_param)]
    d_cds = [0 for i in range(len_seq_param)]
    e_cds = [0 for i in range(len_seq_param)]
    f_cds = [0 for i in range(len_seq_param)]
    g_cds = [0 for i in range(len_seq_param)]
    h_cds = [0 for i in range(len_seq_param)]
    i_cds = [0 for m in range(len_seq_param)]
    j_cds = [0 for i in range(len_seq_param)]
    k_cds = [0 for i in range(len_seq_param)]
    l_cds = [0 for i in range(len_seq_param)]
    m_cds = [0 for i in range(len_seq_param)]
    n_cds = [0 for i in range(len_seq_param)]
    o_cds = [0 for i in range(len_seq_param)]
    p_cds = [0 for i in range(len_seq_param)]
    q_cds = [0 for i in range(len_seq_param)]
    r_cds = [0 for i in range(len_seq_param)]
    s_cds = [0 for i in range(len_seq_param)]
    t_cds = [0 for i in range(len_seq_param)]
    u_cds = [0 for i in range(len_seq_param)]
    v_cds = [0 for i in range(len_seq_param)]
    w_cds = [0 for i in range(len_seq_param)]
    x_cds = [0 for i in range(len_seq_param)]
    y_cds = [0 for i in range(len_seq_param)]
    z_cds = [0 for i in range(len_seq_param)]
    aa_cds = [0 for i in range(len_seq_param)]
    ab_cds = [0 for i in range(len_seq_param)]
    ac_cds = [0 for i in range(len_seq_param)]
    ad_cds = [0 for i in range(len_seq_param)]
    ae_cds = [0 for i in range(len_seq_param)]

    ll = len(parameter_map_tss['normalized_params_map'])

    for mm in range(ll):
        a = np.add(a, parameter_map_tss['normalized_params_map'][mm]['a'])
        b = np.add(b, parameter_map_tss['normalized_params_map'][mm]['b'])
        c = np.add(c, parameter_map_tss['normalized_params_map'][mm]['c'])
        d = np.add(d, parameter_map_tss['normalized_params_map'][mm]['d'])
        e = np.add(e, parameter_map_tss['normalized_params_map'][mm]['e'])
        f = np.add(f, parameter_map_tss['normalized_params_map'][mm]['f'])
        g = np.add(g, parameter_map_tss['normalized_params_map'][mm]['g'])
        h = np.add(h, parameter_map_tss['normalized_params_map'][mm]['h'])
        i = np.add(i, parameter_map_tss['normalized_params_map'][mm]['i'])
        j = np.add(j, parameter_map_tss['normalized_params_map'][mm]['j'])
        k = np.add(k, parameter_map_tss['normalized_params_map'][mm]['k'])
        l = np.add(l, parameter_map_tss['normalized_params_map'][mm]['l'])
        m = np.add(m, parameter_map_tss['normalized_params_map'][mm]['ma'])
        n = np.add(n, parameter_map_tss['normalized_params_map'][mm]['n'])
        o = np.add(o, parameter_map_tss['normalized_params_map'][mm]['o'])
        p = np.add(p, parameter_map_tss['normalized_params_map'][mm]['p'])
        q = np.add(q, parameter_map_tss['normalized_params_map'][mm]['q'])
        r = np.add(r, parameter_map_tss['normalized_params_map'][mm]['r'])
        s = np.add(s, parameter_map_tss['normalized_params_map'][mm]['s'])
        t = np.add(t, parameter_map_tss['normalized_params_map'][mm]['t'])
        u = np.add(u, parameter_map_tss['normalized_params_map'][mm]['u'])
        v = np.add(v, parameter_map_tss['normalized_params_map'][mm]['v'])
        w = np.add(w, parameter_map_tss['normalized_params_map'][mm]['w'])
        x = np.add(x, parameter_map_tss['normalized_params_map'][mm]['x'])
        y = np.add(y, parameter_map_tss['normalized_params_map'][mm]['y'])
        z = np.add(z, parameter_map_tss['normalized_params_map'][mm]['z'])
        aa = np.add(aa, parameter_map_tss['normalized_params_map'][mm]['aa'])
        ab = np.add(ab, parameter_map_tss['normalized_params_map'][mm]['ab'])
        ac = np.add(ac, parameter_map_tss['normalized_params_map'][mm]['ac'])
        ad = np.add(ad, parameter_map_tss['normalized_params_map'][mm]['ad'])
        ae = np.add(ae, parameter_map_tss['normalized_params_map'][mm]['ae'])

    cl=len(parameter_map_cds['normalized_params_map'])

    for mm in range(cl):
        a_cds = np.add(a_cds, parameter_map_cds['normalized_params_map'][mm]['a'])
        b_cds = np.add(b_cds, parameter_map_cds['normalized_params_map'][mm]['b'])
        c_cds = np.add(c_cds, parameter_map_cds['normalized_params_map'][mm]['c'])
        d_cds = np.add(d_cds, parameter_map_cds['normalized_params_map'][mm]['d'])
        e_cds = np.add(e_cds, parameter_map_cds['normalized_params_map'][mm]['e'])
        f_cds = np.add(f_cds, parameter_map_cds['normalized_params_map'][mm]['f'])
        g_cds = np.add(g_cds, parameter_map_cds['normalized_params_map'][mm]['g'])
        h_cds = np.add(h_cds, parameter_map_cds['normalized_params_map'][mm]['h'])
        i_cds = np.add(i_cds, parameter_map_cds['normalized_params_map'][mm]['i'])
        j_cds = np.add(j_cds, parameter_map_cds['normalized_params_map'][mm]['j'])
        k_cds = np.add(k_cds, parameter_map_cds['normalized_params_map'][mm]['k'])
        l_cds = np.add(l_cds, parameter_map_cds['normalized_params_map'][mm]['l'])
        m_cds = np.add(m_cds, parameter_map_cds['normalized_params_map'][mm]['ma'])
        n_cds = np.add(n_cds, parameter_map_cds['normalized_params_map'][mm]['n'])
        o_cds = np.add(o_cds, parameter_map_cds['normalized_params_map'][mm]['o'])
        p_cds = np.add(p_cds, parameter_map_cds['normalized_params_map'][mm]['p'])
        q_cds = np.add(q_cds, parameter_map_cds['normalized_params_map'][mm]['q'])
        r_cds = np.add(r_cds, parameter_map_cds['normalized_params_map'][mm]['r'])
        s_cds = np.add(s_cds, parameter_map_cds['normalized_params_map'][mm]['s'])
        t_cds = np.add(t_cds, parameter_map_cds['normalized_params_map'][mm]['t'])
        u_cds = np.add(u_cds, parameter_map_cds['normalized_params_map'][mm]['u'])
        v_cds = np.add(v_cds, parameter_map_cds['normalized_params_map'][mm]['v'])
        w_cds = np.add(w_cds, parameter_map_cds['normalized_params_map'][mm]['w'])
        x_cds = np.add(x_cds, parameter_map_cds['normalized_params_map'][mm]['x'])
        y_cds = np.add(y_cds, parameter_map_cds['normalized_params_map'][mm]['y'])
        z_cds = np.add(z_cds, parameter_map_cds['normalized_params_map'][mm]['z'])
        aa_cds = np.add(aa_cds, parameter_map_cds['normalized_params_map'][mm]['aa'])
        ab_cds = np.add(ab_cds, parameter_map_cds['normalized_params_map'][mm]['ab'])
        ac_cds = np.add(ac_cds, parameter_map_cds['normalized_params_map'][mm]['ac'])
        ad_cds = np.add(ad_cds, parameter_map_cds['normalized_params_map'][mm]['ad'])
        ae_cds = np.add(ae_cds, parameter_map_cds['normalized_params_map'][mm]['ae'])

    a = [i / ll for i in a]
    b = [i / ll for i in b]
    c = [i / ll for i in c]
    d = [i / ll for i in d]
    e = [i / ll for i in e]
    f = [i / ll for i in f]
    g = [i / ll for i in g]
    h = [i / ll for i in h]
    i = [m / ll for m in i]
    j = [i / ll for i in j]
    k = [i / ll for i in k]
    l = [i / ll for i in l]
    m = [i / ll for i in m]
    n = [i / ll for i in n]
    o = [i / ll for i in o]
    p = [i / ll for i in p]
    q = [i / ll for i in q]
    r = [i / ll for i in r]
    s = [i / ll for i in s]
    t = [i / ll for i in t]
    u = [i / ll for i in u]
    v = [i / ll for i in v]
    w = [i / ll for i in w]
    x = [i / ll for i in x]
    y = [i / ll for i in y]
    z = [i / ll for i in z]
    aa = [i / ll for i in aa]
    ab = [i / ll for i in ab]
    ac = [i / ll for i in ac]
    ad = [i / ll for i in ad]
    ae = [i / ll for i in ae]
    tss_list = [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z, aa, ab, ac, ad, ae]

    a_cds = [i / cl for i in a_cds]
    b_cds = [i / cl for i in b_cds]
    c_cds = [i / cl for i in c_cds]
    d_cds = [i / cl for i in d_cds]
    e_cds = [i / cl for i in e_cds]
    f_cds = [i / cl for i in f_cds]
    g_cds = [i / cl for i in g_cds]
    h_cds = [i / cl for i in h_cds]
    i_cds = [m / cl for m in i_cds]
    j_cds = [i / cl for i in j_cds]
    k_cds = [i / cl for i in k_cds]
    l_cds = [i / cl for i in l_cds]
    m_cds = [i / cl for i in m_cds]
    n_cds = [i / cl for i in n_cds]
    o_cds = [i / cl for i in o_cds]
    p_cds = [i / cl for i in p_cds]
    q_cds = [i / cl for i in q_cds]
    r_cds = [i / cl for i in r_cds]
    s_cds = [i / cl for i in s_cds]
    t_cds = [i / cl for i in t_cds]
    u_cds = [i / cl for i in u_cds]
    v_cds = [i / cl for i in v_cds]
    w_cds = [i / cl for i in w_cds]
    x_cds = [i / cl for i in x_cds]
    y_cds = [i / cl for i in y_cds]
    z_cds = [i / cl for i in z_cds]
    aa_cds = [i / cl for i in aa_cds]
    ab_cds = [i / cl for i in ab_cds]
    ac_cds = [i / cl for i in ac_cds]
    ad_cds = [i / cl for i in ad_cds]
    ae_cds = [i / cl for i in ae_cds]
    cds_list = [a_cds, b_cds, c_cds, d_cds, e_cds, f_cds, g_cds, h_cds, i_cds, j_cds, k_cds, l_cds, m_cds, n_cds, o_cds,
                p_cds, q_cds, r_cds, s_cds, t_cds, u_cds, v_cds, w_cds, x_cds, y_cds, z_cds, aa_cds, ab_cds, ac_cds,
                ad_cds, ae_cds]

    new_plot_dir = os.path.join("plotting/plots/tetra", name)
    try:
        os.mkdir(new_plot_dir)
    except FileExistsError as error:
        print(error)

    for m in range(len(tss_list)):
        plt.clf()
        plt.xlim([0, 1000])
        plt.xticks([i for i in range(0, 1100, 100)], [i for i in range(-500, 600, 100)])
        plt.plot(tss_list[m])
        plt.plot(cds_list[m])
        plt.gcf()
        plt.savefig(str(new_plot_dir)+"/" + str(m) + ".png")

    # print("DONE!!!!!!!")

def plotting_tetra(tetra_tss_seq_map, tetra_cds_seq_map, name):
    print("PLOTTTING!!!!!!!!!!!")

    ll = len(tetra_tss_seq_map)

    l = [0 for i in range(974)]
    m = [0 for i in range(974)]
    n = [0 for i in range(974)]
    o = [0 for i in range(974)]
    p = [0 for i in range(974)]
    q = [0 for i in range(974)]

    l_cds = [0 for i in range(973)]
    m_cds = [0 for i in range(973)]
    n_cds = [0 for i in range(973)]
    o_cds = [0 for i in range(973)]
    p_cds = [0 for i in range(973)]
    q_cds = [0 for i in range(973)]

    for i in range(ll):
        try:
            l = np.add(l,tetra_tss_seq_map[i]['l'])
            m = np.add(m,tetra_tss_seq_map[i]['m'])
            n = np.add(n,tetra_tss_seq_map[i]['n'])
            o = np.add(o,tetra_tss_seq_map[i]['o'])
            p = np.add(p, tetra_tss_seq_map[i]['p'])
            q = np.add(q, tetra_tss_seq_map[i]['q'])
        except :
            print("error at:", i)

    cl=len(tetra_cds_seq_map)

    for i in range(cl):
        l_cds = np.add(l_cds,tetra_cds_seq_map[i]['l'])
        m_cds = np.add(m_cds,tetra_cds_seq_map[i]['m'])
        n_cds = np.add(n_cds,tetra_cds_seq_map[i]['n'])
        o_cds = np.add(o_cds,tetra_cds_seq_map[i]['o'])
        p_cds = np.add(p_cds,tetra_cds_seq_map[i]['p'])
        q_cds = np.add(q_cds,tetra_cds_seq_map[i]['q'])

    l = [i/ll for i in l]
    m = [i / ll for i in m]
    n = [i / ll for i in n]
    o = [i / ll for i in o]
    p = [i / ll for i in p]
    q = [i / ll for i in q]
    tss_list = [l, m, n, o, p, q]

    l_cds = [i / cl for i in l_cds]
    m_cds = [i / cl for i in m_cds]
    n_cds = [i / cl for i in n_cds]
    o_cds = [i / cl for i in o_cds]
    p_cds = [i / cl for i in p_cds]
    q_cds = [i / cl for i in q_cds]
    cds_list = [l_cds, m_cds, n_cds, o_cds, p_cds, q_cds]

    new_plot_dir = os.path.join("../plotting/plots/tetra/", name)
    Path(new_plot_dir).mkdir(parents=True, exist_ok=True)

    for m in range(len(tss_list)):
        plt.clf()
        plt.xlim([0, 1000])
        plt.xticks([i for i in range(0, 1100, 100)], [i for i in range(-500, 600, 100)])
        plt.plot(tss_list[m])
        plt.plot(cds_list[m])
        plt.gcf()
        plt.savefig(str(new_plot_dir)+"/" + str(m) + ".png")

    import sys
    sys.exit()
