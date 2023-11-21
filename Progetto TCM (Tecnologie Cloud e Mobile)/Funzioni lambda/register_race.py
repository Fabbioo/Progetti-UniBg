import boto3
import random

def lambda_handler(event, context):
    
    # Funzione per la registrazione di una gara all'interno di un database DynamoDB.
    # Parametri che devono essere memorizzati:
    # _ id della gara (univoco per ogni gara);
    # _ data della gara;
    # _ nome della gara;
    # _ email del gestore della gara.
    
    # Fase 1: raccolgo i dati dai parametri di input dell'URL.
    race_name = event['queryStringParameters']['race_name']
    race_date = event['queryStringParameters']['race_date'] # Formato: GG-MM-AAAA, ma non è obbligatorio
    email = event['queryStringParameters']['email']
    
    # Fase 2: creo l'id da associare alla gara definita nella fase 1.
    # Questo serve per fare in modo che ogni gara venga identificata da un id.
    race_id = str(random.randint(0, 1000) + random.randint(0, 1000))
        # race_id è la somma di due numeri interi casuali compresi tra 0 e 1000.
        # La somma viene effettuata per ridurre al minimo le probabilità che
        # venga generato una gara con id già presente nel bucket.
    
    # Fase 3: inserisco in una tabella i dati della fase 1 e della fase 2.
    # La tabella in questione si chiama "Palinsesto" e ha come chiave l'id.
    dynamodb = boto3.resource("dynamodb")
    tablePalinsesto = dynamodb.Table("Palinsesto")
    tablePalinsesto.put_item(
        Item = {
            'Id': race_id,
            'Data': race_date,
            'Nome': race_name,
            'Email': email
        }
    )
    
    # Fase 4: creo il token da restituire al gestore della gara.
    token = str(int(random.random() * 1000000))
        # token è un numero compreso tra 0 e 1 moltiplicato per 1000000.
        # Il cast a int serve a eliminare la parte decimale del numero randomico.
        # Il cast a str serve perchè la chiave della tabella è di tipo Stringa.
   
    # Fase 5: inserisco il token nella tabella Autenticazione.
    tableAutenticazione = dynamodb.Table("Autenticazione")
    tableAutenticazione.put_item(
        Item = {
            'Token': token,
            'Nome': race_name,
            'Cognome': race_name,
            'Anno di nascita': 3000 # Valore di default arbitrario
        }
    )
    
    # Fase 6: ritorno al gestore della gara sia il token sia l'id della gara.
    res = {
        'statusCode': 200,
        'body': "token: " + token + "\nrace_id: " + race_id
    }
    
    return res