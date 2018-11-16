#!/bin/bash
if [ -z "$1" ]; then
	echo "Error: invalid sequence ID" >&2
	echo ""
	echo "usage: $0 seqid"
	echo ""
	exit 1
fi

myid=$1
mydb=/usr/users/celldev/luf/database/chr3DL_Bacs_v2.fa
abintio=/usr/users/celldev/luf/temp/Annotation/wheat/1.augustus_training/chr3DLaugustustrainingby3B/run20150709/wheat.augustus.abinitio.20150709.sorted.gff.gz
est=/usr/users/celldev/luf/temp/Annotation/wheat/3.cdna/wheat.cdna.denovo.MIPS3DL.TriFLDB.exonerate.noRM70perc.sorted.gff.gz
ncbi=/usr/users/celldev/luf/temp/Annotation/wheat/3.cdna/ori/ncbi.unigene481773.exonerate.sorted.gff.gz
trinity=/usr/users/celldev/luf/temp/Annotation/wheat/3.cdna/ori/wheat3DL.trinity.exonerate.sorted.gff.gz
mips=/usr/users/celldev/luf/temp/Annotation/wheat/6.MIPSannotation_v1/Ta_3DL-BACs_PGSBv1.0_HighConf_19May2015.sorted.gff3.gz
protein=/usr/users/celldev/luf/temp/Annotation/wheat/4.protein/wheat.protein.gth.sorted.gff.gz
evm=/usr/users/celldev/luf/temp/Annotation/wheat/7.evm/wheat.evm.20150714.NoEmpty.sorted.gff.gz

echo "### 0. Cleaning..."
rm ./*.gz ./*.tbi ./1.*

echo "### 1. fasta..."
fasta.extractor.sh -i $mydb -s $myid -o 1.$myid.fasta
if [ -s "1.$myid.fasta" ]; then
	samtools faidx 1.$myid.fasta
else
	echo "Error: Invalid ID, empty fasta..." >&2
	exit 1
fi

echo "### 2. AB initio..."
tabix $abintio $myid | bgzip > 2.$myid.abinitio.gff.gz
if [ -s "2.$myid.abinitio.gff.gz" ]; then
	tabix -p gff 2.$myid.abinitio.gff.gz
fi

echo "### 3. EST..."
tabix $est $myid | bgzip > 3.$myid.est.gff.gz
if [ -s "3.$myid.est.gff.gz" ]; then
	tabix -p gff 3.$myid.est.gff.gz
fi

echo "### 4. NCBI..."
tabix $ncbi $myid | bgzip > 4.$myid.ncbi.gff.gz
if [ -s "4.$myid.ncbi.gff.gz" ]; then
	tabix -p gff 4.$myid.ncbi.gff.gz
fi

echo "### 5. Trinity..."
tabix $trinity $myid | bgzip > 5.$myid.trinity.gff.gz
if [ -s "5.$myid.trinity.gff.gz" ]; then
	tabix -p gff 5.$myid.trinity.gff.gz
fi

echo "### 6. MIPS..."
tabix $mips $myid | bgzip > 6.$myid.mips.gff.gz
if [ -s "6.$myid.mips.gff.gz" ]; then
	tabix -p gff 6.$myid.mips.gff.gz
fi

echo "### 7. Protein..."
tabix $protein $myid | bgzip > 7.$myid.protein.gff.gz
if [ -s "7.$myid.protein.gff.gz" ]; then
	tabix -p gff 7.$myid.protein.gff.gz
fi

echo "### 8. EVM..."
tabix $evm $myid | bgzip > 8.$myid.evm.gff.gz
if [ -s "8.$myid.evm.gff.gz" ]; then
	tabix -p gff 8.$myid.evm.gff.gz
fi

#art 1.$myid.fasta +2.$myid.abintio.gff.gz +3.$myid.est.gff.gz +4.$myid.ncbi.gff.gz +5.$myid.trinity.gff.gz +6.$myid.mips.gff.gz +7.$myid.protein.gff.gz +8.$myid.evm.gff.gz


exit 0
