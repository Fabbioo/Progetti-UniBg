import json
import boto3

def lambda_handler(event, context):
   
    # Funzione che visualizza tutte le gare disponibili memorizzati in un database DynamoDB
   
    # Fase 1: accedo alla tabella Palinsesto.
    dynamodb = boto3.resource("dynamodb")
    tablePalinsesto = dynamodb.Table("Palinsesto")
    
    # Fase 2: creo un dizionario contenente tutti i dati presenti nel palinsesto.
    response = tablePalinsesto.scan() # response è di tipo dict.
    
    # Fase 3: creo una lista in cui posizionare tutti i dizionari presenti in response['Items'].
    list = []
    for index, evento in enumerate(response['Items']):
        # response['Items'] è una lista, non un dizionario come nel caso di response.
        # evento è un dizionario.
        list.append(evento) # Appendo i dizionari in coda ad una lista.
    
    return {
        'statusCode': 200,
        'body': (json.dumps(list, indent = 4)) # l'output è di tipo str
    }