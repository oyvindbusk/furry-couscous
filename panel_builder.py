import argparse
import csv
import collections
import httplib2 as http
import shutil
import json
import time
import os
import sys
import pandas as pd
import urllib
pd.options.mode.chained_assignment = None
try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse
    from urllib.parse import quote


class Main_method:
    '''
    This class takes a file with gene names as input:
    "
    GJB1
    PMP22
    "
x = df.groupby("name2")["chrom"].apply(','.join).reset_index()
x[x.chrom.str.join(' ').str.contains('M')]
    For each gene name it checks against genenames.org to see if it is valid. If it is not, it wont continue.

    When all genes are valid, a file will be saved that can be used as import to genetikkportalene.no.

    It will also produce bedfiles containing intervals for all exons in the entire gene bp on each exon. The bed containing exon intervals will not emit any genes with no coding exons, but they will be present in the full_gene.bed-file.


    Refseq select: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20201022/GCF_000001405.25_GRCh37.p13/


    '''
    def __init__(self):
        self.genelist       = args.genelist
        self.name           = args.name
        self.genes          = self.read_genes()
        self.genecount      = len(self.genes)
        self.duplicates     = self.find_duplicates()
        self.refseq         = "downloads/RefseqSelect_GRCh37.txt"
        self.refseqAll      = "downloads/RefseqAll_GRCh37.txt"
        self.genCode        = "downloads/GenCodeBasicV36lift37.txt"
        self.notfound       = []
        self.genedf         = pd.read_csv("downloads/hgnc_complete_set.txt", delimiter="\t", usecols=[0, 1], low_memory=False)
        self.dtypes         = {'#bin':int,'name':str,'chrom':str,'strand':str,'txStart':int,'txEnd':int,'cdsStart':int,'cdsEnd':int,'exonCount':int,'exonStarts':str,'exonEnds':str, 'score':int,   'name2':str ,  'cdsStartStat':str ,  'cdsEndStat':str ,'exonFrames':str}

    def run_local_query(self, gene):
        # Load protein-coding_gene.txt into a df
        result = self.genedf.loc[self.genedf['symbol'] == gene]
        if len(result) != 0:
            return (gene, result['hgnc_id'].to_string(index=False).lstrip())
        else:
            print "!! %s were not found in the database." % (gene)
            print "Checking previous symbols:"
            self.run_query(gene)
            self.notfound.append(gene)

    def limit_query(self):
        ''' Runs query trough this because of server load limit'''
        geneIDList = []
        for count, gene in enumerate(self.genes):
            geneIDList.append(g.run_local_query(gene))
        if len(self.notfound) != 0:
            print "The following genes were not found in the database: %s" % (','.join(self.notfound))
        return [gene for gene in geneIDList if gene]

    def read_genes(self):
        genes = []
        with open(self.genelist, 'rb') as g:
            reader = csv.reader(g)
            for row in reader:
                if row:
                    genes.append(row[0])
        return genes

    def run_query(self, gene):
        ''' If a gene symbol is not in the downloaded list, check if it is in previous symbol in genenames.org'''
        headers = {'Accept': 'application/json',}
        uri = 'http://rest.genenames.org'
        path = '/search/prev_symbol/' + gene
        target = urlparse(uri+path)
        method = 'GET'
        body = ''
        h = http.Http()
        response, content = h.request(target.geturl(),method,body,headers)
        notfound = []
        if response['status'] == '200':
        # parse content with the json module 
            data = json.loads(content)
            if len(data['response']['docs']) != 0:
                print "\t * %s was found in previous symbol. Could the new symbol be:%s" % (gene, data['response']['docs'][0]['symbol'])
                return data['response']['docs'][0]['symbol']
            else:
                notfound.append(gene)
        else:
            print 'Error detected: ' + response['status']
        
    def run_query_prev_symbol(self, gene):
        ''' Find the previous symbol of a gene name'''
        headers = {'Accept': 'application/json',}
        uri = 'http://rest.genenames.org'
        path = '/fetch/symbol/' + gene
        target = urlparse(uri+path)
        method = 'GET'
        body = ''
        h = http.Http()
        response, content = h.request(target.geturl(),method,body,headers)
        notfound = []
        if response['status'] == '200':
        # parse content with the json module 
            data = json.loads(content)
            if len(data['response']['docs']) != 0:
                if 'prev_symbol' in data['response']['docs'][0].keys():
                    return data['response']['docs'][0]['prev_symbol']
                else:
                    print "\t\t\t * No previous symbol found"          
        else:
            print 'Error detected: ' + response['status']
        

    def find_duplicates(self):
        ''' Find duplicates in gene list: '''
        return [item for item, count in collections.Counter(self.genes).items() if count > 1]

    def find_duplicates_in_bed(self, df):
        return df[df.duplicated(subset=['name2'])]    


    def explode_df(self, df):
        '''
        Takes the exonstarts and exonends fields, and splits them up one for each row:
        '''
        listOfTuples = []
        for index, row in df[df["exonCount"] > 1 ].iterrows(): #skipping genes with only 1 exon
            if row['strand'] == "+":
                for c, (i,j) in enumerate(zip(row['exonStarts'].split(',')[:-1], row['exonEnds'].split(',')[:-1])):
                    listOfTuples.append(tuple(row[0:2]) + (row[2], ) + tuple(row[3:8]) + (row[8], i , j) + (row[11], row[12]+'_E' + str(c + 1)) + tuple(row[13:15])+ (row[15].split(',')[c],))
            else:
                exoncount = row[8]
                for c, (i,j) in enumerate(zip(row['exonStarts'].split(',')[:-1], row['exonEnds'].split(',')[:-1])):
                    listOfTuples.append(tuple(row[0:2]) + (row[2], ) + tuple(row[3:8]) + (row[8], i , j) + (row[11], row[12]+'_E' + str(exoncount - (c))) + tuple(row[13:15])+ (row[15].split(',')[c],))


        # Adding genes with one exon:
        df[["exonStarts", "exonEnds","exonFrames"]] = df.loc[:,["exonStarts", "exonEnds","exonFrames"]].replace(",$","",regex=True) # remove trailing comma
        listOfTuples.append(df[df["exonCount"] == 1].to_records(index=False)[0])

        result_df = pd.DataFrame(listOfTuples, columns = df.columns)
        # Remove those exons with exonframes -1 (complete UTR) unless the entire gene is UTR (e.g. SNORD118, TERC)
        result_df = result_df[(result_df.exonFrames != "-1") | (result_df.exonCount == 1)]

        # Remove parts of the exons that are not part of CDS
        # make a mask of all rows where the exonstarts are not within interval cdsstart-cdsend
        startmask = ((~result_df.exonStarts.astype(int).between(result_df.cdsStart,result_df.cdsEnd)) & ((result_df.exonCount.astype(int) != 1) & (result_df.exonFrames != "-1")))
        endmask   = ((~result_df.exonEnds.astype(int).between(result_df.cdsStart,result_df.cdsEnd)) & ((result_df.exonCount.astype(int) != 1) & (result_df.exonFrames != "-1")))
        result_df.loc[startmask, 'exonStarts'] = result_df.cdsStart # For those, set exonstarts = cddsstarts
        result_df.loc[endmask, 'exonEnds']     = result_df.cdsEnd # For those, set exonstarts = cddsstarts
        return result_df

    def save_as_genebed(self, df):
        ''' Formats df into bed format one line for each gene'''
        df = df[['chrom','txStart','txEnd','name','score','strand']]
        df.to_csv("%s/%s_full_gene.BED" % (self.name,self.name), sep="\t", header=False, index=False)

    def save_as_exonbed(self, df):
        ''' Formats df into bed format one line pr exon'''
        df = df[['chrom','exonStarts','exonEnds','name','name2','score','strand']]
        df.to_csv("%s/%s.BED" % (self.name,self.name), sep="\t", header=False, index=False)

    def find_noncoding(self,df):
        ''' find genes that are entirely noncoding e.g. SNORD118'''
        dfy = df[~df.duplicated('name2')]
        return dfy[(dfy["cdsStart"] - dfy["cdsEnd"] == 0)]["name2"]

    def find_mito(self, df):
        dfi = df.groupby("name2")["chrom"].apply(','.join).reset_index()
        return dfi[dfi.chrom.str.join(' ').str.contains('M')]      
    
    def noncoding(self, df):
        ''' If exoncount is 1 and frame is -1 (e.g. TERC, SNORD118) use cdsstarts/ends instead of exonstart/end in bed'''



    def make_bed(self):
        ''' Take the output file from the previous step, get the HGNC-IDs and make a bed out of it'''
        df          = pd.read_csv(self.refseq,      delimiter="\t", low_memory=False, dtype=self.dtypes)
        df_all      = pd.read_csv(self.refseqAll,   delimiter="\t", low_memory=False, dtype=self.dtypes)
        df_genCode  = pd.read_csv(self.genCode,     delimiter="\t", low_memory=False, dtype=self.dtypes)
        all_genes   = pd.DataFrame()
        old_genes   = []  # List of genes where previous symbol was used instead of the updated one - must be checked manually
        notfound    = []  # List of genes not found at all with previsous or updated symbol in refseq select or all
        for g in self.genes:
            if len(df.loc[df['name2'] == g]) != 0:
                all_genes = pd.concat([all_genes, df.loc[df['name2'] == g]])
                print "* Success: %s" % (g)
            else:
                print "Gene %s is not found in Refseq Select, Trying refseq all..." % (g)
                if len(df_all.loc[df_all['name2'] == g]) != 0:
                    print "\t * Success: %s" % (g)
                    # Since not using the refseqSelect file, must check that only one transcript pr gene:
                    all_genes = pd.concat([all_genes, df_all.loc[df_all['name2'] == g]])
                else:
                    print "Gene %s is not found in Refseq All, Trying Gencode..." % (g)
                    if len(df_genCode.loc[df_genCode['name2'] == g]) != 0:
                        print "\t * Success: %s" % (g)
                        all_genes = pd.concat([all_genes, df_genCode.loc[df_genCode['name2'] == g]])
                    else:
                        prevgene = ""
                        print "\t\t\t * Trying with previous symbol instead" # Some of the gene identifiers in the refseq-files are not updated. Therefore get the previous symbol from genenames.org. These must be ckecked manually
                        prevgene = self.run_query_prev_symbol(g) # in case multiple previous names
                        if not prevgene:
                            print "\t\t\t\t * Not found using previous - Abort - Abort"
                            sys.exit()
                        for c, pg in enumerate(prevgene):
                            if len(df_all.loc[df_all['name2'] == pg]) != 0:
                                all_genes = pd.concat([all_genes, df_all.loc[df_all['name2'] == pg]])
                                old_genes.append((g, pg))
                                print "\t\t\t\t * Success: %s" % (g)
                                break
                            else:
                                if c - len(prevgene) == 1:
                                    print "\t\t\t\t * Not found"
                                    notfound.append(g)
                                else:
                                    print "\t\t\t\t * Not found, trying next of the previous symbols"

        # Sort dataframe
        all_genes.sort_values(by=['chrom','exonStarts','exonEnds'], inplace=True)

        # Warn check where old symbols are used:
        if len(old_genes) != 0:
            print "######################################################################################################"
            print "Some genes were not found in Refseq at all. They were found when using old symbols. Check these manually !!!"
            print "The genes are: %s" % ([ "%s: %s" % (x,y) for x,y in zip([i for sub in old_genes for i in sub][0::2], [i for sub in old_genes for i in sub][1::2]) ])
        
        # Warn if some genes not found at all:
        if len(notfound) != 0:
            print "######################################################################################################"
            print "Some genes were not found at all !!!"
            print "The genes are: %s" % (", ".join(notfound))

        # Warn if duplicates are present in output:
        # Duplicated gene names:
        if len(self.find_duplicates_in_bed(all_genes)) != 0:
            print "######################################################################################################"
            print "Some genes have duplicate entries. This applies to the following:"
            dupl = self.find_duplicates_in_bed(all_genes)["name2"]
            dupllist = all_genes[all_genes["name2"].isin( list(dupl))][['name','name2', 'exonCount', 'chrom','exonStarts', 'exonEnds','txStart','txEnd']]
            dupllist["size"] = abs(dupllist["txStart"] - dupllist["txEnd"])
            print dupllist
            # Write to file:
            all_genes[all_genes["name2"].isin( list(dupl))][['name2', 'chrom','exonStarts', 'exonEnds']].to_csv("%s/%s_duplicated.regions.txt" % (self.name,self.name), sep="\t", header=True, index=False)
       
        ## Print genes that are entirely non-conding e.g SNORD118
        if len(self.find_noncoding(all_genes)) != 0:
            print "######################################################################################################"
            print "The following genes are entirely nonconding and needs to looked at manually:"
            for i in self.find_noncoding(all_genes):
                print i
        if len(self.find_mito(all_genes)):
            print "######################################################################################################"
            print "The following looks like mitochondrial genes?"
            print(self.find_mito(all_genes).to_string()) 

        print "######################################################################################################"

        # Save pr gene bed
        self.save_as_genebed(all_genes)
        # Save pr exon bed
        self.save_as_exonbed(self.explode_df(all_genes))


if __name__ =="__main__":

    # Argument parsing
    parser = argparse.ArgumentParser(description='Takes a file with gene names as inputm and queries each name against genenames.org, to check if the genename is present', epilog='')
    parser.add_argument('-g','--genelist', help='Input file', required=True)
    parser.add_argument('-n','--name', help='Name of panel', required=True)
    args = parser.parse_args()

    # Init main
    g = Main_method()

    # Create folders
    os.makedirs(g.name)

    # find Duplicates:
    if len(g.duplicates) != 0:
        print "######################################################################################################"
        print "Duplicates are present in list: %s" % (','.join(g.duplicates))
        print "Consider removing them before moving on."

    # Number of genes
    print "######################################################################################################"
    print "Panel name: %s" % (args.name)
    print "Making input files for genetikkportalen.no"
    print "######################################################################################################"
    print "There are %s genes in the genefile." % (g.genecount)
    formatted_genes = g.limit_query()
    print "There are %s genes in the output file" % (len(formatted_genes))
    print "Done"
    
    # Write to file genetikkportalen output
    with open("%s/%s.genetikkportalimport.txt" % (g.name, g.name),'wb') as out:
        csv_out=csv.writer(out, delimiter=";")
        csv_out.writerow(["HGNC_symbol","HGNC_id"])
        csv_out.writerows(formatted_genes)

    # Rename genelist-file and put in outputfolder
    shutil.copy(g.genelist, "%s/%s.txt" % (g.name, g.name))
    


    # Make bedfiles
    print "######################################################################################################"
    print "Making bedfiles:"
    print "######################################################################################################"
    g.make_bed()
    print "Done"
    print "######################################################################################################"




'''
Check the following genes:

!! ARMC4 were not found in the database.
Checking previous symbols:
         * ARMC4 was found in previous symbol. Could the new symbol be:ODAD2
!! C12orf65 were not found in the database.
Checking previous symbols:
         * C12orf65 was found in previous symbol. Could the new symbol be:MTRFR
!! CCDC114 were not found in the database.
Checking previous symbols:
         * CCDC114 was found in previous symbol. Could the new symbol be:ODAD1
!! CCDC151 were not found in the database.
Checking previous symbols:
         * CCDC151 was found in previous symbol. Could the new symbol be:ODAD3
!! LRRC6 were not found in the database.
Checking previous symbols:
         * LRRC6 was found in previous symbol. Could the new symbol be:DNAAF11
!! TTC25 were not found in the database.
Checking previous symbols:
         * TTC25 was found in previous symbol. Could the new symbol be:ODAD4
!! WDR34 were not found in the database.
Checking previous symbols:
         * WDR34 was found in previous symbol. Could the new symbol be:DYNC2I2




'''


# Run like:
