def empirical_hpd(values, conf=0.05):
    """
    Assuming a **unimodal** distribution, returns the 0.95 highest posterior
    density (HPD) interval for a set of samples from a posterior distribution.
    Adapted from ``emp.hpd`` in the "TeachingDemos" R package (Copyright Greg
    Snow; licensed under the Artistic License).
    """
#    conf = min([conf, 1.0 - conf])

    conf = 1.0 - conf
    n = len(values)
    nn = int(round(n * conf))
    x = sorted(values)
    xx = []
    if nn == 0:
        raise ValueError("Sample size too small: %s" % len(values))
    for i in range(nn):
        Z1 = x[n-(nn-i)]
        Z2 = x[i]
        #print "==>", Z1, Z2
        xx.append(Z1 - Z2)
    m = min(xx)
    nnn = xx.index(m)
    #print "n=", n
    #print "nn=", nn
    #print "xx=", xx
    #print "m=", m
    #print "nnn=", n
    return [x[nnn], x[n-nn+nnn]]

def empirical_mode(values,frac_1sigma=5):

  # first work out what precision level we care about - what is the 1 sigma width?


  sigma1 = empirical_hpd(values,0.6827)

  careabout = np.diff(sigma1)/frac_1sigma

  start = 0.6827

  mode = np.mean(sigma1)

#  print "let's go!",sigma1

  for i in range(0,100):
    start = start/1.1
    new_sigma1 = empirical_hpd(values,start)
    new_mode = np.mean(new_sigma1)
#    print 'new_sigma',new_sigma1
#    print new_mode
    if abs(abs(mode-new_mode)) < careabout:
      if (new_mode > sigma1[0]) & (new_mode < sigma1[1]):
        return new_mode
      else:
#        print 'something wrong with mode estimation, using median'
        return new_mode
#        return np.median(values)
    mode = new_mode

  if (new_mode > sigma1[0]) & (new_mode < sigma1[1]):
    return new_mode
  else:
    print('something wrong with mode estimation, using median')
    return np.median(values)
