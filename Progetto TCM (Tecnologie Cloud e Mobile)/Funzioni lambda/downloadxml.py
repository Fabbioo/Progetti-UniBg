import boto3
import base64

def lambda_handler(event, context):
   
   # Funzione per il download di un file xml dal bucket
    
    file_name = event['queryStringParameters']['id']
    s3 = boto3.client('s3')
    response = s3.get_object(
        Bucket = "myawsbuckettcm",
        Key = file_name + ".xml"
    )
    
    image_to_be_downloaded = response['Body'].read()
    
    return {
        'statusCode': 200,
        'body': base64.b64encode(image_to_be_downloaded),
        'isBase64Encoded': True
    }
    