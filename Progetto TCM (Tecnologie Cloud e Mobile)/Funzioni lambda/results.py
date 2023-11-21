import json
import xml.etree.ElementTree as ET
import boto3
import datetime

def lambda_handler(event, context):
    
    # Funzione per la visualizzazione dei risultati di una gara.
    
    # Parametri di input:
    # _ id della gara;
    # _ categoria della gara di cui estrarre i risultati.
    
    # Parametri di output:
    # _ posizione;
    # _ id dell'atleta;
    # _ nome dell'atleta;
    # _ cognome dell'atleta;
    # _ data e ora di inizio della gara;
    # _ data e ora di fine della gara;
    # _ tempo impiegato in minuti e secondi;
    # _ nome del team dell'atleta.
    
    race_id_param = event["queryStringParameters"]["id"]
    race_id_param = race_id_param + ".xml"
    categoria_param = event["queryStringParameters"]["class"]

    BUCKET_NAME = "myawsbuckettcm"
    s3 = boto3.resource('s3')

    s3obj = s3.Object(BUCKET_NAME, race_id_param)
    
    xmlstring = s3obj.get()['Body'].read().decode('ISO-8859-1')
    
    root = ET.fromstring(xmlstring)
    
    _ns = {
        "":'http://www.orienteering.org/datastandard/3.0'
    }
        
    class_names_elements = root.findall("./ClassResult/Class[Name = '" + str(categoria_param) + "']/../PersonResult", _ns)
        
    lista = []
    posizione = 0
    for y in class_names_elements:
        posizione += 1
        if y.find('./Person//Id', _ns).text == None:
            id = "id non specificato"
        else:
            id = y.find('./Person//Id', _ns).text
        name = y.find('./Person//Name/Given', _ns).text
        surname = y.find('./Person//Name/Family', _ns).text
        startTime = y.find('.//Result/StartTime', _ns).text
        finishTime = y.find('.//Result/FinishTime', _ns).text
        
        if y.find('.//Result/Time', _ns).text == None:
            time = "sconosciuto"
        else:
            time = y.find('.//Result/Time', _ns).text
            time = str(datetime.timedelta(seconds=int(time))) # Conversione del tempo da secondi a minuti
        organiz = y.find('.//Organisation/Name', _ns).text
        results = {}
        results["posizione"] = str(posizione)
        results["id"] = id
        results["atleta"] = name + " " + surname
        results["tempoInizio"] = startTime
        results["tempoFine"] = finishTime
        results["tempo"] = time
        results["organizzazione"] = organiz
        lista.append(results)
    return {
        'statusCode': 200,
        'body': json.dumps(lista)
    }
