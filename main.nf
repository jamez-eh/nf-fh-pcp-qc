nextflow.preview.dsl=2


process fastqc {
	container 'nfcore/rnaseq:1.4.2'

	input:
	tuple val(sampleID), val(type), val(kitID), val(patientID), path(fastq1), path(fastq2)

	output:
    	file "*_fastqc.{zip,html}"


	"""
	    fastqc --quiet --threads 3 ${fastq1} ${fastq2}

	"""
}

process fastq_pair {
        container 'chgyi/fastq_pair'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 50

        input:
        tuple val(sampleID), val(type), val(kitID), val(patientID), file(fastq1), file(fastq2)

        output:
        tuple file("${sampleID}_R1.fq.paired.fq.gz"), file("${sampleID}_R2.fq.paired.fq.gz")


        """
	gunzip -c ${fastq1} > ${sampleID}_R1.fq
	gunzip -c ${fastq2} > ${sampleID}_R2.fq
	fastq_pair -t 100000 ${sampleID}_R1.fq ${sampleID}_R2.fq
	gzip *.fq
        """
}

workflow fastqc_wf {

         take: 
	       fastqs
         main:
	       fastqc(fastqs)
	 emit:
		qc = fastqc.out		
}

workflow fix_fq_wf {
	  take:
		fastqs
	
	  main:
		fastq_pair(fastqs)
	  emit:
		fixed = fastq_pair.out
}


	



workflow {
	 
        fq_ch = Channel
            .fromPath(params.input_csv)
            .splitCsv(header:true)
            .map{ row-> tuple(row.sampleID, row.type, row.kitID, row.patientID, file(row.R1), file(row.R2)) }
	

	main:
		//fastqc_wf(fq_ch)
		fix_fq_wf(fq_ch)
	
	publish: 
		fix_fq_wf.out.fixed to : "/fh./fixed/"   
		//fastqc_wf.out.qc to : "./qc/"
		
}
