process  getSra {
  input:
  publishDir "/home/orfeas/Documents/sraQuery/output/"
    path csv_file
  output:
    env ids 

  script:
  """
  awk -F "\"*,\"*" '{print \$2}' $csv_file
  """
}

process  fetch {
  input:
  publishDir "/home/orfeas/Documents/sraQuery/output/"
    val id_str

  script:
  """
  prefetch $id_str.id
  """
}