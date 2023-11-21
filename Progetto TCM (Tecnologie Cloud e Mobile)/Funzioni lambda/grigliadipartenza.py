import json
import xml.etree.ElementTree as ET
import boto3

def lambda_handler(event, context):

    # Funzione per estrapolare dai file xml la griglia di partenza degli atleti.
    
    # Parametri di input:
    # _ id della gara;
    # _ categoria della gara di cui estrarre la griglia.
    
    # Parametri di output:
    # _ una lista contenente i seguenti dati:
    #   _ id dell'atleta;
    #   _ nome dell'atleta;
    #   _ cognome dell'atleta;
    #   _ data e ora di partenza.

    s3 = boto3.resource('s3')
    try:
        raceid = event['queryStringParameters']['id']
        category = event['queryStringParameters']['class']
        bucket_name = 'myawsbuckettcm'
        object_name = str(raceid) + '.xml'
        obj = s3.Object(bucket_name, object_name)
        xmlstr = obj.get()['Body'].read()
    except s3.meta.client.exceptions.NoSuchKey:
        return {
            'statusCode': 404,
            'body': json.dumps("Errore: id non trovato.")
        }
    except KeyError:
        return {
            'statusCode': 400,
            'body': json.dumps("Errore: id non specificato.")
        }

    # STEP 2: parsing dell'XML.
    root = ET.fromstring(xmlstr)
    _ns = {'': 'http://www.orienteering.org/datastandard/3.0'}
    class_names_elements = root.findall("./ClassResult/Class[Name='" + str(category) + "']/../PersonResult", _ns)
    
    lista = []
    # Query xPath per l'estrazione delle informazioni necessarie al completamento della lista
    for y in class_names_elements:
        if y.find('./Person//Id', _ns).text == None:
            id = "id non specificato"
        else:
            id = y.find('./Person//Id', _ns).text
        name = y.find('./Person//Name/Given', _ns).text
        surname = y.find('./Person//Name/Family', _ns).text
        starttime = y.find('.//Result/StartTime',_ns).text
        results = {}
        results["id"] = id
        results["name"] = name
        results["surname"] = surname
        results["starttime"] = starttime
        lista.append(results)
    return {
        'statusCode': 200,
        'body': json.dumps(lista)
    }