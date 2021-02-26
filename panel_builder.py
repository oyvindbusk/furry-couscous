import argparse
import csv
import collections
import httplib2 as http
import json
import time
import os
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
        self.genes = self.read_genes()
        self.genecount = len(self.genes)
        self.duplicates = self.find_duplicates()
        
    
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

            #print 'Symbol:' + data['response']['docs'][0]['symbol']
        else:
            print 'Error detected: ' + response['status']

        if len(notfound) != 0:
            print "Some genes were not found in the database: %s" % (','.join(notfound))













    def limit_query(self):
        ''' Runs query trough this because of server load limit'''
        geneIDList = []
        for count, gene in enumerate(self.genes):
            geneIDList.append(g.run_query(gene))
            if count % 3 == 0:
                time.sleep(1)
        
        return [gene for gene in geneIDList if gene] 

    def get_bed(self, gene):
        ''' query bed regions for each gene    ''' 
        headers = {'Accept': 'application/json',}
        uri = 'http://mygene.info'
        path = '/v3/gene/' + gene
        method = 'GET'
        target = urlparse(uri+path)
        body = ''
        h = http.Http()
        response, content = h.request(target.geturl(),method,body,headers)
        if response['status'] == '200':
        # parse content with the json module 
            data = json.loads(content)
            print json.dumps(data['exons_hg19'], indent=4, sort_keys=True)


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

    def make_downloads(self):
        ''' Make folder and download files needed for exons '''
        try:
            os.makedirs('downloads')
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        if path.exists("downloads")

if __name__ =="__main__":
    parser = argparse.ArgumentParser(description='Takes a file with gene names as inputm and queries each name against genenames.org, to check if the genename is present', epilog='')
    parser.add_argument('-g','--genelist', help='Input file', required=True)
    #parser.add_argument('-v', '--vcf', help='Input vcf from sequencing', required=True)
    args = parser.parse_args()
    g = Main_method()
    
    # Duplicates:
    if len(g.duplicates) != 0:
        print "Duplicates are present in list: %s" % (','.join(g.duplicates))
        print "Consider removing them before moving on."
    # Stats:
    print "There are %s genes in the genefile." % (g.genecount)
    print "\n"
    #g.run_query("GJB1111")
    #g.get_bed("1017")
    print g.limit_query()
    
    # else:
    #     output_data = converter.intersect()
    #     print('Intersects the data with Bedtools\n')
    #     print('.......')
    #     converter.save_file(output_data)
    #     print('Converts output files to the correct input format for annovar.')
    #     print('.......')
    #     converter.make_annovar_file()
    #     print('Minor adjustment and combination of files.')
    #     print('.......')
    #     converter.refine_annovar_file()
    #     converter.run_annovar()
    #     converter.cleanup()
