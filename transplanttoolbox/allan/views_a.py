from django.shortcuts import render
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.schemas import SchemaGenerator
from rest_framework.permissions import AllowAny
from rest_framework import status, generics
from rest_framework_swagger.renderers import SwaggerUIRenderer, OpenAPIRenderer
import hla
import conversion_functions
from conversion_functions import convert_allele_to_ag
from hla import allele_truncate
from . import serializers
from rest_framework_swagger.views import get_swagger_view



class AlleleApiView(generics.GenericAPIView):
    """Returns UNOS antigen for an allele. Enter IPD-IMGT/HLA allele fully resolved or upto second field with expression characters. The results yield UNOS antigen equivalency and Bw4/6 epitope. """
    serializer_class = serializers.AlleleSerializer
    permission_classes = [AllowAny,]
        #def get(self, response, format=None):
        #obj = ["UNOS antigen mapping for WHO HLA alleles"]
        #return Response({'Web Services': obj})

    def post(self, request, format=None):
        """Returns UNOS antigen for an allele.Enter IPD-IMGT/HLA allele fully resolved or upto second field with expression characters. The results yield UNOS antigen equivalency and Bw4/6 epitope. """
        """parameters:
            allele:string
            """
        serializer = serializers.AlleleSerializer(data=request.data)    

        if serializer.is_valid(raise_exception=True):
            allele = serializer.data.get('allele')
            ag = conversion_functions.convert_allele_to_ag(allele)
            result = {'Antigen': ag}
            return Response(result, status=status.HTTP_200_OK)
        else:
            return Response({"Error": "Check if allele is IMGT/HLA"}, status=status.HTTP_400_BAD_REQUEST)
                #serializer.errors, status=status.HTTP_400_BAD_REQUEST)    
