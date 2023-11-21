import boto3
from boto3.dynamodb.conditions import Key

def lambda_handler(event, context):
    
    # Fase 1: verifico che il token sia valido.
    dynamodb = boto3.resource('dynamodb')
    table = dynamodb.Table("Autenticazione")
    id = event['queryStringParameters']['token']
    response = table.query(KeyConditionExpression = Key('Token').eq(id))
    Count = int(response['Count'])
    
    # Fase 2: se il token è valido procedo con l'upload nel bucket del file xml.
    if Count != 0:   
        if id == response['Items'][0]['Token']:
            xmlstring = event['body']
            object_name = event['queryStringParameters']['raceid'] + ".xml"
            s3 = boto3.resource('s3')
            obj = s3.Object("myawsbuckettcm", object_name)
            obj.put(Body = xmlstring)
            res = {
                'statusCode': 200,
                'body': 'File caricato correttamente.'
            }
    else:
        res = {
            'statusCode': 200,
            'body': 'Utente non autorizzato al caricamento.'
        }

    return res