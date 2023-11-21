import json
import xml.etree.ElementTree as ET
import boto3

def lambda_handler(event, context):

    # Funzione per estrapolare dai file xml il nome di tutte le gare disponibili
    
    # STEP 1: recupero della strina XML nel bucket S3 corrispondente all'ID passato come parametro
    s3 = boto3.resource('s3')
    try:
        raceid = event['queryStringParameters']['id']
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
    class_names_elements = root.findall('./ClassResult/Class/Name', _ns)
    classes = [x.text for x in class_names_elements]

    return {
        'statusCode': 200,
        'body': json.dumps(classes)
    }