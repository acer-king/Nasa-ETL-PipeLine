import ast
import boto3
import botocore
import os
import logging
from botocore.exceptions import ClientError
from . import re_snowflake

logging.basicConfig(
    format="%(asctime)s  %(levelname)s: %(message)s", level=logging.INFO
)


def get_secret(secret_name: str, region_name: str = "us-east-1") -> dict:
    """
    Utility function to retrieve secrets from AWS secret manager as a dictionary
    :param secret_name: AWS secret name
    :param region_name: AWS reqion name
    :return: Dictionary containing secrets
    """

    client = boto3.client(service_name="secretsmanager", region_name=region_name)

    # In this sample we only handle the specific exceptions for the 'GetSecretValue' API.
    # See https://docs.aws.amazon.com/secretsmanager/latest/apireference/API_GetSecretValue.html
    # We rethrow the exception by default.

    try:
        get_secret_value_response = client.get_secret_value(SecretId=secret_name)
    except ClientError as e:
        if e.response["Error"]["Code"] == "DecryptionFailureException":
            # Secrets Manager can't decrypt the protected secret text using the provided KMS key.
            # Deal with the exception here, and/or rethrow at your discretion.
            raise e
        elif e.response["Error"]["Code"] == "InternalServiceErrorException":
            # An error occurred on the server side.
            # Deal with the exception here, and/or rethrow at your discretion.
            raise e
        elif e.response["Error"]["Code"] == "InvalidParameterException":
            # You provided an invalid value for a parameter.
            # Deal with the exception here, and/or rethrow at your discretion.
            raise e
        elif e.response["Error"]["Code"] == "InvalidRequestException":
            # You provided a parameter value that is not valid for the current state of the resource.
            # Deal with the exception here, and/or rethrow at your discretion.
            raise e
        elif e.response["Error"]["Code"] == "ResourceNotFoundException":
            # We can't find the resource that you asked for.
            # Deal with the exception here, and/or rethrow at your discretion.
            raise e
    else:
        # Decrypts secret using the associated KMS CMK.
        # Depending on whether the secret is a string or binary, one of these fields will be populated.
        if "SecretString" in get_secret_value_response:
            secret: dict = ast.literal_eval(get_secret_value_response["SecretString"])
            return secret
        else:
            return None
    finally:
        del client


def get_credentials(secrets: dict = {}) -> dict:
    """Get Credentials

    Get credentials by implementing the following rules of retrieval:
    1. If parameter secrets is provided, is a dict and not empty, return as-is.
    2. Otherwise, try getting credentials from environment variables.
    3. If that's still unsucessful, try with the AWS secrets manager via get_secret().

    For step (2), this function supports the following environment variables:
    - Snowflake: SF_USER, SF_PASS and SF_HOST (standardized in output as SNOWFLAKE_USER,
    SNOWFLAKE_PASS and SNOWFLAKE_ACCOUNT, respectively)
    - Quandl: QUANDL_API_KEY

    :param secrets: Dictionary containing secrets. Empty dictionary {} by default.

    :return: Credentials dictionary. If no data can be found, an empty dictionary.
    """

    # If param secrets is a non-empty dictionary, return as is
    if isinstance(secrets, dict) and len(secrets):
        logging.info("Returning secrets parameter as-is.")

        return secrets

    # Initialize helper vars
    credentials = {}
    template_error_msg = "Error in getting credentials: \n {}"

    # Try to retrieve from environment variables
    logging.info("Getting credentials from environment variables.")
    try:
        # NOTE: Add here any credential we need. At the moment we are supporting
        # Snowflake connection settings and the API key for Quandl
        credentials = {
            "SNOWFLAKE_USER": os.environ["SF_USER"],
            "SNOWFLAKE_PASS": os.environ["SF_PASS"],
            "SNOWFLAKE_ACCOUNT": os.environ["SF_HOST"].replace(
                ".snowflakecomputing.com", ""
            ),
            "QUANDL_API_KEY": os.environ["QUANDL_API_KEY"],
        }

        # Return here if retrieval was successful
        return credentials
    except KeyError as e:
        logging.info("Credentials not found in environment variables: \n {}".format(e))
    except Exception as e:
        logging.exception(template_error_msg.format(e))

    # If impossible to get from env vars, try the AWS secrets manager
    logging.info("Getting credentials from AWS secrets manager.")
    try:
        credentials = get_secret(secret_name="redas/prod", region_name="us-east-1")
    except Exception as e:
        logging.exception(template_error_msg.format(e))

    if not secrets:
        logging.warning("No credentials found, returning an empty dictionary.")

    return credentials


def s3_object_exists(bucket: str, key: str) -> bool:
    """
    Check is an s3 object exists.
    :param bucket: String
                    Bucket name
    :param key: String
                    s3 key
    :return: Boolean
                True: if s3 object exists,
                False: if s3 object doesn't exist.
    """
    s3 = boto3.client('s3')
    obj_exists = True
    try:
        s3.head_object(Bucket=bucket, Key=key)
    except botocore.exceptions.ClientError:
        obj_exists = False
    finally:
        del s3
    return obj_exists