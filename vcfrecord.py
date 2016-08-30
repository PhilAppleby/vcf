# Helper methods for parsing VCF records - very basic just now
class VCFrecord():
  def __init__(self):
    self.first_genotype_idx = 9
    self.chr_idx = 0
    self.posn_idx = 1
    self.rsid_idx = 2
    self.ref_idx = 3
    self.alt_idx = 4
    self.fmt_idx = 8

  def get_data_array(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    return (full_record.split('\t'))

  def get_chr(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    return (full_record.split('\t')[self.chr_idx])

  def get_chr_from_array(self, data):
    """
    The argument is a split VCF record
    """
    return (data[self.chr_idx])

  def get_var_id(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    return (full_record.split('\t')[self.rsid_idx])

  def get_var_id_from_array(self, data):
    """
    The argument is a split VCF record
    """
    return (data[self.rsid_idx])

  def get_alleles(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    data = full_record.split('\t')
    return (data[self.ref_idx], data[self.alt_idx])

  def get_alleles_from_array(self, data):
    """
    The argument is a split VCF record
    """
    return (data[self.ref_idx], data[self.alt_idx])

  def get_posn(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    return (full_record.split('\t')[self.posn_idx])

  def get_posn_from_array(self, data):
    """
    The argument is a split VCF record
    """
    return (data[self.posn_idx])

  def get_posn_as_int(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    return (int(full_record.split('\t')[self.posn_idx]))

  def get_posn_from_array_as_int(self, data):
    """
    The argument is a split VCF record
    """
    return (int(data[self.posn_idx]))

  def get_fmts(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    data = full_record.split('\t')
    return (data[self.fmt_idx].split(':'))

  def get_fmts_from_array(self, data):
    """
    The argument is a split VCF record
    """
    return (data[self.fmt_idx].split(':'))

  def get_fmt_indices(self, fmts, req_fmts):
    """
    Arguments are 2 lists of strings
    """
    idxs = {}
    for idx, fmt in fmts:
      if fmt in req_fmts:
        idxs[fmt] = idx
    return idxs

  def set_fmts_in_array(self, data, fmts):
    """
    The argument are a split VCF record and an array of formats
    """
    data[self.fmt_idx] = ":".join(fmts)
    return data

  def get_prfx_sfx(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    rec = full_record.split('\t')
    return (rec[:self.first_genotype_idx], rec[self.first_genotype_idx:])

  def get_prfx_sfx_from_array(self, data):
    """
    The argument is an split VCF record
    """
    return (data[:self.first_genotype_idx], data[self.first_genotype_idx:])

  def get_call_rate(self, full_record):
    """
    """
    prfx,sfx = self.get_prfx_sfx(full_record)
    num_genotypes = len(sfx)
    call_count = 0
    for geno in sfx:
      try:
        (varcall, probs) = geno.split(":")
      except:
        continue
      if varcall != "./.":
        call_count += 1

    return (float(call_count) / num_genotypes)

  def get_call_rate_from_array(self, data):
    """
    """
    prfx,sfx = self.get_prfx_sfx_from_array(data)
    num_genotypes = len(sfx)
    call_count = 0
    for geno in sfx:
      try:
        (varcall, probs) = geno.split(":")
      except:
        continue
      if varcall != "./.":
        call_count += 1

    return (float(call_count) / num_genotypes)

  def get_genotype_counts(self, full_record):
    """
    """
    geno_counts = {}
    called_samp_count = 0
    samp_count = 0

    prfx, sfx = self.get_prfx_sfx(full_record)

    for geno in sfx:
      (varcall, probs) = geno.split(":")
      geno_counts[varcall] = geno_counts.get(varcall, 0) + 1
      samp_count += 1
      if varcall != "./.":
        called_samp_count += 1

    return geno_counts, called_samp_count, samp_count

  def get_genotype_counts_from_array(self, data):
    """
    """
    geno_counts = {}
    called_samp_count = 0
    samp_count = 0

    prfx, sfx = self.get_prfx_sfx_from_array(data)

    for geno in sfx:
      (varcall, probs) = geno.split(":")
      geno_counts[varcall] = geno_counts.get(varcall, 0) + 1
      samp_count += 1
      if varcall != "./.":
        called_samp_count += 1

    return geno_counts, called_samp_count, samp_count

  def get_allele_counts_from_array(self, data):
    """
    """
    homref_count = 0
    het_count = 0
    homalt_count = 0
    nc_count = 0
    miss_count = 0 # can be "." or ""

    prfx, sfx = self.get_prfx_sfx_from_array(data)

    for geno in sfx:
      if geno == "." or geno == "":
        miss_count += 1
      else:
        (varcall, probs) = geno.split(":")
        if varcall == "0/0":
          homref_count += 1
        elif varcall == "1/1":
          homalt_count += 1
        elif varcall == "./.":
          nc_count += 1
        else:
          het_count += 1

    return homref_count, het_count, homalt_count, nc_count, miss_count
