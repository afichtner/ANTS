import xml.etree.ElementTree as et

def read_xml(filename):
    
    def recursive_dict(element):
        return element.tag, \
            dict(map(recursive_dict, element)) or element.text
    
    
    doc = et.parse(filename)
    
    return recursive_dict(doc.getroot())



def find_coord(path_to_xml):
    tree=et.parse(path_to_xml)
    root=tree.getroot()
    
    sta=path_to_xml.split('/')[-1].split('.')[1]

    lat=root.find('*//Station').find('Latitude').text
    lon=root.find('*//Station').find('Longitude').text
    return sta, float(lat),float(lon)
    
    
    
    