import re
import uuid

def make_fasta(repo, date):
    fasta = {}
    fasta['datecreated'] = date
    fasta['features'] = repo['features']
    fasta['nsa'] = repo['nsa']
    fasta['nucseq'] = repo['nucseq']
    fasta['parts'] = repo['parts']
    fasta['persons'] = {} # TODO make subjects dict - repo['subjects']
    fasta['repo'] = repo
    return fasta


def save_part(fasta, partname, partauthorname, partdescription, partsequence, datecreated):
    partsequencenonl = re.sub(r'[^AGCTagct]+', '', partsequence)
    partid = uuid.uuid4()

    seqid = partid
    nucseq = [nuc for nuc in fasta['nucseq'] if nuc['idnucseq'] == seqid]
    if not nucseq:
        nucseq = {}
        nucseq['idnucseq'] = seqid
        nucseq['sequence'] = partsequencenonl
        fasta['nucseq'].append(nucseq)

    part = {}
    part['name'] = partname
    part['authorid'] = partauthorname
    part['description'] = partdescription
    part['idpart'] = partid
    part['isbasic'] = 0
    part['riskgroup'] = 0
    part['nucseq'] = nucseq
    part['datecreated'] = datecreated
    fasta['parts'].append(part)

    feature  = {}
    feature['name'] = partname
    feature['authorid'] = partauthorname
    feature['idfeature'] = partid
    feature['riskgroup'] = 0
    feature['nucseq'] = nucseq
    feature['datecreated'] = datecreated
    fasta['features'].append(feature)

    return part