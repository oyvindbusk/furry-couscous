import argparse
import csv
import collections
import httplib2 as http
import json
import time
import os
import pandas as pd
pd.options.mode.chained_assignment = None
try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse









class Main_method:
    ''' 
    This class takes a file with gene names as input:
    "
    GJB1
    PMP22
    "
    
    For each gene name it checks against genenames.org to see if it is valid. If it is not, it wont continue. 

    When all genes are valid, a file will be saved that can be used as import to genetikkportalene.no.

    It will also produce bedfiles containing intervals for all exons in the entire gene +- 10 bp on each exon
   
    Refseq select: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20201022/GCF_000001405.25_GRCh37.p13/


    '''
    def __init__(self):
        self.genelist = args.genelist
        self.name = args.name
        self.genes = self.read_genes()
        self.genecount = len(self.genes)
        self.duplicates = self.find_duplicates()
        self.refseq = "downloads/RefseqSelect_GRCh37.txt"
        self.notfound = []
    
    def run_query(self, gene):
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
                print "%s;%s" % (gene, data['response']['docs'][0]['hgnc_id'])
                return (gene, data['response']['docs'][0]['hgnc_id'])
            else:
                notfound.append(gene)
                self.notfound.append(gene)
        else:
            print 'Error detected: ' + response['status']
        if len(notfound) != 0:
            print "!! %s were not found in the database." % (','.join(notfound))

    def limit_query(self):
        ''' Runs query trough this because of server load limit'''
        geneIDList = []
        for count, gene in enumerate(self.genes):
            geneIDList.append(g.run_query(gene))
            if count % 3 == 0:
                time.sleep(1)
        if len(self.notfound) != 0:
            print "The following genes were not found in the database: %s" % (','.join(self.notfound))
        return [gene for gene in geneIDList if gene] 

    def read_genes(self):
        genes = []
        with open(self.genelist, 'rb') as g:
            reader = csv.reader(g)
            for row in reader:
                genes.append(row[0])
        return genes

    def find_duplicates(self):
        ''' Find duplicates in gene list: '''
        return [item for item, count in collections.Counter(self.genes).items() if count > 1]

    def explode_df(self, df):
        ''' 
        Takes the exonstarts and exonends fields, and splits them up one for each row:
        '''
        
        listOfTuples = []
        for index, row in df.iterrows():
            for c, (i,j) in enumerate(zip(row['exonStarts'].split(',')[:-1], row['exonEnds'].split(',')[:-1])):
                listOfTuples.append(tuple(row[0:2]) + (row[2].replace('chr', ''), ) + tuple(row[3:8]) + (1, i , j) + (row[11],row[12]+'_E' + str(c + 1)) + tuple(row[13:15])+ (row[15].split(',')[c],))
        return pd.DataFrame(listOfTuples, columns = df.columns)

    def format_to_bed(self, df):
        ''' Formats df into bed format and adds 10 pb on each interval'''
        return df[['chrom','exonStarts','exonEnds','name','score','strand']]

    def save_as_bed(self, df):
        
        df.to_csv("%s/%s.BED" % (self.name,self.name), sep="\t", header=False, index=False)    
    
    def make_bed(self):
        df = pd.read_csv(self.refseq, delimiter="\t")
        all_genes = pd.DataFrame()
        for g in self.genes:
            all_genes = pd.concat([all_genes, df.loc[df['name2'] == g]])
            if len(df.loc[df['name2'] == g]) == 0:
                print "Gene %s is not found in Refseq Select, consider adding manually" % (g)
        self.save_as_bed(self.format_to_bed(self.explode_df(all_genes)))


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
        print "Duplicates are present in list: %s" % (','.join(g.duplicates))
        print "Consider removing them before moving on."
    
    # Number of genes
    print "There are %s genes in the genefile." % (g.genecount)
    formatted_genes = g.limit_query()
    print "There are %s genes in the output file" % (len(formatted_genes))
    # Write to file
    with open("%s/%s.genetikkportalimport.txt" % (g.name, g.name),'wb') as out:
        csv_out=csv.writer(out, delimiter=";")
        csv_out.writerow(["HGNC_symbol","HGNC_id"])
        csv_out.writerows(formatted_genes)
    # Make bedfiles
    print "Making bedfiles:"
    g.make_bed()
    print "Done"







# Run like: 
# python panel_builder.py -g testfiles/genepanels/Schwan_v2_210121.txt -n Schwan_v2_210121