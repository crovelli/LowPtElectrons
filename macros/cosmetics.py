class RangesByName(object):
   def __init__(self, rlist):
      self._rlist_ = rlist #list of (ending, range)

   def get(self, val, default=None):
      for ending, vrange in self._rlist_:
         if val.endswith(ending):
            return vrange
      return default

   def __getitem__(self, val):
      self.get(val)

ranges = RangesByName([
   ('_cluster_deta', (-2, 2)),
   ('_pt', (0, 15)),
   ('_eta' , (-3, 3)),
   ('_inp' , (0, 20)),
   ('_outp' , (0, 10)),
   ('_chi2red' , (0, 6)),
   ('_Deta' , (0, 0.2)),
   ('_Dphi' , (-0.2, 0.2)),
   ('_nhits' , (0, 50)),
   ('_p' , (0, 20)),
   ('_cluster_e', (0, 20)),
   ('_cluster_ecorr', (0, 20)),
   ('_cluster_eta', (-3, 3)),
   ('_cluster_deta', (-1.5, 1.5)),
   ('_cluster_dphi', (-1.5, 1.5)),
   ## ('_cluster_covEtaEta', (-0.5, 0.5)),
   ## ('_cluster_covEtaPhi', (-0.5, 0.5)),
   ## ('_cluster_covPhiPhi', (-0.5, 0.5)),
])

beauty = {
    'gen_pt' : r'p$_T$(gen)', 
    'gen_eta' : r'$\eta$(gen)', 
    'dxy_err' : r'$\sigma$(dxy)', 
    'nhits' : r'\# of hits',
    'trk_pt' : r'p$_T$(ktf track)', 
    'trk_eta' : r'$\eta$(ktf track)', 
    'log_trkpt' : r'log$_{10}$(p$_T$)(ktf track)', 
    'trk_inp' : r'p$_{in}$(ktf track)',
    'trk_outp' : r'p$_{out}$(ktf track)', 
    'trk_eta' : r'$\eta$(ktf track)', 
    'trk_ecal_Deta': '$\Delta\eta$(ECAL, ktf track)',
    'trk_ecal_Dphi' : '$\Delta\varphi$(ECAL, ktf track)',
    'e_over_p' : 'E/p', 
    'trk_chi2red' : '$\chi^2$(ktf track)/ndf', 
    'gsf_dpt' : r'p$_T$(gsf track)',
    'trk_gsf_chiratio' : '$\chi^2$(gsf track)/$\chi^2$(ktf track)', 
    'gsf_chi2red' : '$\chi^2$(gsf track)/ndf', 
    'xy_sig' : r'$\sigma$(dxy)/dxy',
}
