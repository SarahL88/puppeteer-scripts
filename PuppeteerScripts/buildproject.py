import json

import repository

def generate_build_request(repo, concentration, concentration_unit, volume_unit, authorid):
    partssofar = {}
    vectorssofar = {}
    partslibrary = []
    vectorlibrary = []
    designlist = []

    # TODO original was: for part in user-selected parts:
    for gff3design in repo['gff3designs']:

        constituents = repository.get_constituent_parts(repo, gff3design['name'])

        partnames = {}
        position = 0
        for part in constituents:
            if part['name'] not in partssofar:
                pp = make_part(repo, part, concentration, concentration_unit, volume_unit)
                partslibrary.append(pp)
                partssofar[part['name']] = part
            partnames[position] = part['name'] + '-' + str(part['idpart'])
            position += 1

        design = {}
        design['name'] = gff3design['name']
        design['partPositionMap'] = partnames
        designlist.append(design)

        # TODO - original processed vectors in 'repository.getvectorsbypart(repo, part)'
        # But the backbone vector (first vector in vectors.csv) is not associated with a part.
        # So, here, we use the list of all vectors (repo['vectors']).
        allvectors = repo['vectors']
        for v in allvectors:
            if v['name'] not in vectorssofar:
                vectorssofar[v['name']] = v
                vec = make_vector(repo, v, concentration, concentration_unit, volume_unit)
                vectorlibrary.append(vec)

        backbone_vector = allvectors[0]
        design['vectorName'] = backbone_vector['name'] + '-' + str(backbone_vector['idvector'])

        buildrequest = {}
        buildrequest['partSamples'] = partslibrary
        buildrequest['vectorSamples'] = vectorlibrary
        buildrequest['designs'] = designlist
        buildrequest['doBuildabilityVerification'] = False
        buildrequest['lineSeparator'] = '\n'

        parameters = {}
        parameters['volume'] = 20.0
        parameters['buildmethod'] = 'Modular Cloning v1.0'
        buildrequest['parameters'] = parameters
        buildrequest['ownerUuid'] = authorid

    request = json.dumps(buildrequest, indent=2)
    return request


def make_part(repo, p, concentration, concentration_unit, volume_unit):
    part = {}
    part['name'] = p['name'] + '-' + str(p['idpart'])
    part['sequence'] = p['nucseq']['sequence'].upper()
    part['concentration'] = concentration
    part['concentrationUnit'] = concentration_unit
    part['volumeUnit'] = volume_unit
    nsa = repository.get_annotations_by_family(repo, p['nucseq'], 'overhang')

    l = []
    r = []
    if nsa[0]['startx'] > nsa[1]['startx']:
        l.append(nsa[1]['startx'])
        l.append(nsa[1]['endx'])
        r.append(nsa[0]['startx'])
        r.append(nsa[0]['endx'])
    else:
        l.append(nsa[0]['startx'])
        l.append(nsa[0]['endx'])
        r.append(nsa[1]['startx'])
        r.append(nsa[1]['endx'])

    overhang = {}
    overhang['fivePrimeEnd'] = l[0]
    overhang['fivePrimeStart'] = l[1]
    overhang['threePrimeEnd'] = r[0]
    overhang['threePrimeStart'] = r[1]
    part['overhangs'] = overhang
    return part


def make_vector(repo, v, concentration, concentration_unit, volume_unit):
    vector = {}
    vector['name'] = v['name'] + '-' + str(v['idvector'])
    vector['sequence'] = v['nucseq']['sequence'].upper()
    vector['concentration'] = concentration
    vector['concentrationUnit'] = concentration_unit
    vector['volumeUnit'] = volume_unit

    locs = repository.get_moclo_vector_digestion_locations(repo, v);
    endofvectorprefix = locs[1]
    vectorprefix = v['nucseq']['sequence'][0:endofvectorprefix + 1]
    startofvectorprefix = locs[2]
    vectorsuffix = v['nucseq']['sequence'][startofvectorprefix:]

    overhang = {}
    overhang['fivePrimeEnd'] = locs[0]
    overhang['fivePrimeStart'] = locs[1]
    overhang['threePrimeEnd'] = locs[1] + 1
    overhang['threePrimeStart'] = locs[1] + 1 + (locs[3] - locs[2])

    vector['sequence'] = (vectorprefix + vectorsuffix).upper()
    vector['overhangs'] = overhang
    return vector


def print_nested_dict(d, indent=0):
    # For debugging
    for key, value in d.items():
        if isinstance(value, str):
            print('\t' * indent + str(key) + '\t\t' + str(value))
            continue;
        else:
            print('\t' * indent + str(key))
        if isinstance(value, list):
            for v in value:
                if (isinstance(v, dict)):
                    print_nested_dict(v, indent + 1)
        if isinstance(value, dict):
            print_nested_dict(value, indent + 1)
        else:
            if isinstance(value, str):
                print('\t' * (indent+1) + str(value))
