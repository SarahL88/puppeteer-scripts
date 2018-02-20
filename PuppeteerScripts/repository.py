import uuid
import random
import re

def create_feature(repo, featurename, featuresequence, familyname, date):
    feature = {}
    featureid = uuid.uuid4()
    feature['idfeature']= featureid
    feature['datecreated']=date
    feature['forcolor']=round(1 + random.random() * 13)
    feature['iscds']= ''
    feature['name']= featurename

    nucseq = {}
    nucseq['datecreated'] =date
    nucseq['idnucseq']=featureid
    nucseq['sequence']=featuresequence

    feature['nucseq'] = nucseq
    repo['features'].append(feature)

    family = get_family_by_name(repo, familyname)

    ffxref = {}
    ffxref['family'] = family
    ffxref['feature'] = feature
    ffxref['datecreated'] = date
    ffxref['lastmodified'] = date

    ffxrefpk = {}
    ffxrefpk['featureid'] = featureid;
    ffxrefpk['familyid'] = family['idfamily']
    ffxref['ffxrefpk'] = ffxrefpk

    repo['ffxref'].append(ffxref)

    return feature;



def get_family_by_name(repo, familyname):
    if repo['families']:
        for fam in repo['families']:
            if 'name' in fam:
                if fam['name'] == familyname:
                    return fam;
                else:
                    if fam['name'] == 'User-defined':
                        return fam

    family = {}
    family['idfamily'] = familyname #TODO what should id be
    family['name'] = familyname
    repo['families'].append(family)
    return family



def add_object_to_collection(repo, collectionid, objectid, objecttype, authorid, date):
    cxref = {}
    cxrefpk = {}
    cxrefpk['collectionid']=collectionid
    cxrefpk['objectid']= objectid
    cxref['cxrefpk'] = cxrefpk

    cxref['objecttype'] = objecttype.upper()
    cxref['datecreated'] = date
    cxref['authorid'] = authorid
    cxref['lastmodified'] = date
    cxref['xrefid'] =uuid.uuid4()
    repo['cxref'].append(cxref)


def get_features_by_family_name(repo, familyname):
    family = get_family_by_name(repo, familyname);
    featurefamilies = [ffx for ffx in repo['ffxref'] if ffx['family']['idfamily'] == family['idfamily']]
    features = [f['feature'] for f in featurefamilies]
    return features



def add_feature_to_nucseq(repo, name, nucseq, feature, position, authorid, date):
    nucseqannotation = {}
    nucseqannotation['authorid'] = authorid
    nucseqannotation['datecreated'] = date
    nucseqannotation['feature'] = feature
    nucseqannotation['forwardcolor'] = feature['forcolor']
    nucseqannotation['name'] = feature['name']
    nucseqannotation['nucseq'] = nucseq
    nucseqannotation['reversecolor'] = feature['forcolor']

    if position > len(nucseq['sequence']) - 1 or position < 0:
        raise ValueError('Invalid position ', position,
                         ' provided: Feature ', feature['name'],
                         'referred to in ', name, ' does not occur in its sequence.')

    nucseqannotation['startx'] = position
    nucseqannotation['endx'] = position + len(feature['nucseq']['sequence']) - 1
    nucseqannotation['annotationid'] = uuid.uuid4()

    repo['nsa'].append(nucseqannotation)


def get_annotations_by_family(repo, nucseq, familyname):
    overhanganno = []
    allnucseqnsa = [nsa for nsa in repo['nsa'] if nsa['nucseq']['idnucseq'] == nucseq['idnucseq']]
    for nsa in allnucseqnsa:
        feature = nsa['feature']
        ffx = [ffx for ffx in repo['ffxref'] if ffx['feature']['idfeature'] == feature['idfeature']][0]
        fam = ffx['family']
        if 'name' in fam:
            if fam['name'].lower() == familyname.lower():
                overhanganno.append(nsa)

    return overhanganno



def get_overhang_position_in_plasmid(nucseq, overhangsequence, dir):
    '''
    Gets first occurrence of overhang sequence in the nucseq.
    The overhang sequence is abutted by:
        1) a BsaI-N- site upstream, or
        2) a N-reverse-BsaI downstream.
    '''

    overhangsequence = overhangsequence.strip().lower()

    forwardBsaISitePosStrand = "ggtctc"
    forwardBbsISitePosStrand = "gaagac"
    reverseBsaISitePosStrand = "gagacc"
    forwardBsaIDistance = 1
    forwardBbsIDistance = 2

    maxupstreamlen = max(len(forwardBsaISitePosStrand) + forwardBsaIDistance, len(forwardBbsISitePosStrand) + forwardBbsIDistance)

    circularizedseq = nucseq[len(nucseq) - 1 - (maxupstreamlen - 1):] + nucseq
    circularizedseq = circularizedseq.strip().lower()

    if dir == 'FIVE_PRIME':
        regx = (re.escape(forwardBsaISitePosStrand) + r'[agct]' + re.escape(overhangsequence))
    elif dir == 'THREE_PRIME':
        regx = (re.escape(overhangsequence) + r'[agct]' + re.escape(reverseBsaISitePosStrand))
    else:
        raise ValueError('Invalid DNA Direction: ' + dir)

    match = re.search(regx, circularizedseq, re.IGNORECASE)
    if match:
        return match.start()

    return -1



def get_overhang_position_in_vector(nucseq, overhangsequence):
    '''
    Gets first occurrence of overhang sequence in the nucseq.
    The overhang sequence is abutted by:
        1) a BsaI-N- upstream and NN-BbsI downstream, or
        2) a BsbI-NN-upstream and N-BsaI downstream.
    '''

    overhangsequence = overhangsequence.strip()

    forwardBsaISitePosStrand = "ggtctc"
    forwardBbsISitePosStrand = "gaagac"
    reverseBbsISitePosStrand = "gtcttc"
    reverseBsaISitePosStrand = "gagacc"
    forwardBsaIDistance = 1
    forwardBbsIDistance = 2

    # Get max length of potential abutting strands
    maxupstreamlen = max(len(forwardBsaISitePosStrand) + forwardBsaIDistance,
                         len(forwardBbsISitePosStrand) + forwardBbsIDistance)

    # Concat onto nucseq a substring of nucseq - the substring length is the max abutting strand size
    circularizedseq = nucseq[(len(nucseq) - 1 - maxupstreamlen - 1):] + nucseq;
    circularizedseq = circularizedseq.lower()

    regex = []
    regex.append(re.escape(forwardBbsISitePosStrand) + r'[agct][agct]' + re.escape(overhangsequence) + r'[agct]' + re.escape(reverseBsaISitePosStrand))
    regex.append(re.escape(forwardBsaISitePosStrand) + r'[agct]' + re.escape(overhangsequence) + r'[agct][agct]' + re.escape(reverseBbsISitePosStrand))

    for rx in regex:
        match = re.search(rx, circularizedseq, re.IGNORECASE)
        if match:
            return match.start()

    return -1


def persist_part(repo, partname, partsequence, description, isbasic, authorid, date):
    part = {}
    part['authorid'] = authorid
    part['datecreated'] = date
    part['lastmodified'] = date
    if 'e0030_cd' in partname.lower():      #TODO - Move this to input validation check
        partname = 'Part-E0030m_CD'
    part['name'] = partname
    partid = uuid.uuid4()
    part['idpart'] = partid
    if isbasic:
        part['basic'] = 1
    else:
        part['basic'] = 0
    part['description'] = description

    format = [f for f in repo['formats'] if f['idformat'] == 'edu-bu-synbiotools-format-moclo']
    if len(format) > 0:
        part['format'] = format[0]

    nucseq = {}
    nucseq['datecreated'] = date
    nucseq['lastmodififed'] = date
    nucseq['sequence'] = partsequence
    nucseq['idnucseq'] = partid
    repo['nucseq'].append(nucseq)

    part['nucseq'] = nucseq
    repo['parts'].append(part)

    return part

def persist_plasmid(repo, name, part, vector, authorid, date):
    plasmid = {}
    plasmid['authorid'] = authorid
    plasmid['datecreated'] = date
    plasmid['lastmodified'] = date

    format = [f for f in repo['formats'] if f['idformat'] == 'edu-bu-synbiotools-format-moclo']
    if len(format) > 0:
        plasmid['format'] = format[0]

    plasmid['part'] = part
    plasmid['vector'] = vector
    plasmid['name'] = name
    plasmidid = uuid.uuid4()
    plasmid['idplasmid'] = plasmidid
    repo['plasmids'].append(plasmid)

    # TODO: note - orig Java Puppeteer persists Simple Rich Sequence

    return plasmid


def get_constituent_parts(repo, designname):
    return [cx['parentpart'] for cx in repo['compositexrefs'] if cx['childpart']['name'] == designname]


def get_vectors_by_part(repo, part):
    plasmids = [p for p in repo['plasmids'] if p['part']['idpart'] == part['idpart']]
    vectors = [plasmid['vector'] for plasmid in plasmids]
    return vectors


def get_moclo_vector_digestion_locations(repo, vector):
    nsa = get_annotations_by_family(repo, vector['nucseq'], 'overhang')
    if len(nsa) != 2:
        raise ValueError('Vector ' + vector['name'] + ' contains an unexpected number (' + \
                         len(nsa) + ' of overhang + annotations.')

    if nsa[0]['startx'] < nsa[1]['startx']:
        vectorminanno = nsa[0]
        vectormaxanno = nsa[1]
    else:
        vectorminanno = nsa[1]
        vectormaxanno = nsa[0]

    locations= [vectorminanno['startx'],
                vectorminanno['endx'],
               vectormaxanno['startx'],
                vectormaxanno['endx']]
    return locations


def get_moclo_overhang_annotation(repo, nucseq, dir):
    nsa = get_annotations_by_family(repo, nucseq, 'overhang')
    if len(nsa) < 2:
        raise ValueError('Unexpectedly few annotations were found on this sequence;'
                         'cannot find all overhangs.')
    min = float('inf')
    max = float('-inf')

    for n in nsa:
        if n['startx'] < min:
            min = n['startx']
            minanno = n
        if n['startx'] > max:
            max = n['startx']
            maxanno = n

    if dir == 'FIVE_PRIME':
        return minanno
    elif dir == 'THREE_PRIME':
        return maxanno
    else:
        raise ValueError('Unexpected parameter value in query.')


def find_vectors_by_overhang_resistance(repo, fiveprimeoverhang, threeprimeoverhang, compatibleresistances):
    matching = []
    for vector in repo['vectors']:
        nt = vector['nucseq']
        nsa = get_annotations_by_family(repo, nt, 'resistance')

        resistancefound = False
        for n in nsa:
            if n['feature'] in compatibleresistances:
                resistancefound = True
                break

        overhangsfound = False
        nsa = get_annotations_by_family(repo, nt, 'overhang')
        if nsa[0]['startx'] < nsa[1]['startx']:
            if nsa[0]['feature']['idfeature'] == fiveprimeoverhang['idfeature'] and \
                            nsa[1]['feature']['idfeature'] == threeprimeoverhang['idfeature']:
                        overhangsfound = True
        else:
            if nsa[0]['feature']['idfeature'] == threeprimeoverhang['idfeature'] and \
                            nsa[1]['feature']['idfeature'] == fiveprimeoverhang['idfeature']:
                        overhangsfound = True

        if resistancefound and overhangsfound:
            matching.append(vector)

    return matching


# Functions in original Java Puppeteer repository.java but not here:
# removeallfiltersexceptsubject
# getplasmidsbypart
# getpartsincollection

