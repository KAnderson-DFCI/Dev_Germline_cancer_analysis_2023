version development

workflow VCFfixup {
  input {
    File vcf
    File vcf_idx
    String? pref
  }

  call vcflib_vcffixup {
    input:
    vcf = vcf,
    vcf_idx = vcf_idx,
    pref = pref
  }

  output {
    File fixed_vcf = vcflib_vcffixup.fixed_vcf
    File fixed_vcf_idx = vcflib_vcffixup.fixed_vcf_idx
  }
}

task vcflib_vcffixup {
  input {
      File vcf
      File vcf_idx
      String? pref

      String outFormat = 'z'

      Int memGB = 16
      Int cpu = 6
      Int preemptible = 3
      String vcflib_docker = "kylera/vcflib:latest"
  }

  Int diskGB = ceil(2.5 * size(vcf, 'GB')) + 10

  String suf = if outFormat == "v" then ".vcf" else (
               if outFormat == "z" then ".vcf.gz" else (
               if outFormatf == "u" then ".bcf" else (
               if outFormat == "b" then ".bcf.gz" else ".vcf")))

  String here_vcf = basename(vcf)
  String base = sub(here_bcf, "\\.[bv]cf(\\.gz)?", "")
  String out = if defined(pref) then "~{pref}.~{suf}" else "~{base}.AC_fixed.~{suf}"

  command <<<
    ln -s ~{vcf} ~{here_vcf}
    ln -s ~{vcf_idx} ~{basename(vcf_idx)}

    vcffixup ~{here_vcf} \
    | bcftools view \
      -O ~{outFormat} \
      -o ~{out} \
      --threads 3

    bcftools index -f ~{out} --threads 3
  >>>

  output {
    File fixed_vcf = "~{out}"
    File fixed_vcf_idx = "~{out}.csi"
  }

  runtime {
    disks: "local-disk ~{diskGB} HDD"
    memory: "~{memGB}GB"
    cpu: cpu
    preemptible: preemptible
    docker: vcflib_docker
  }
}