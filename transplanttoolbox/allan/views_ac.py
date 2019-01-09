from django.shortcuts import render
from rest_framework.schemas import SchemaGenerator
from rest_framework.permissions import AllowAny
from rest_framework_swagger.renderers import SwaggerUIRenderer, OpenAPIRenderer
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status, generics
import hla
import conversion_functions
from conversion_functions import convert_allele_to_ag, convert_allele_list_to_ags, gl_string_ags, allele_code_ags
from hla import allele_truncate
from . import serializers
from rest_framework_swagger.views import get_swagger_view


class AlleleCodesApiView(generics.GenericAPIView):
    """Returns antigens for a list of allele codes from 3 through 6 loci (A, B, C, DRB1, DQB1, DRB3/4/5). Enter acronym for the population typed for and the results show the most probable antigens and presences of Bw4 and Bw6 epitopes is indicated. Antigen probabilities for all the possible antigen genotypes mapped to the allele codes are also displayed."""
    serializer_class = serializers.MACSerializer
        #def get(self, response, format=None):
        #obj = ["UNOS antigen mapping for WHO HLA alleles"]
        #return Response({'Web Services': obj})

    def post(self, request):
        """Returns UNOS antigen for a list of  allele codes from 3 through 6 loci (A, B, C, DRB1, DQB1, DRB3/4/5). Enter acronym for the population typed for and the results show the most probable antigens and presences of Bw4 and Bw6 epitopes is indicated. Antigen probabilities for all the possible antigen genotypes mapped to the allele codes are also displayed."""
        serializer = serializers.MACSerializer(data=request.data)    

        if serializer.is_valid():
            allele_codes_list = serializer.data.get('allele_codes_list')
            allele_codes_split = allele_codes_list.split(" ")
            population = serializer.data.get('pop')
            output = conversion_functions.allele_code_ags(allele_codes_split, population)
            print(output)
            ags = output[0::3]
            bw4_6_list = output[1::3]
            probs = output[2::3]
            return Response({'Allele Codes List': allele_codes_list, 'Race': population, 'Antigens': ags, "Bw4/6": bw4_6_list, "Antigen Probabilities": probs})
        else:
            return Response({"Error": "Check if allele is IMGT/HLA"})
                #serializer.errors, status=status.HTTP_400_BAD_REQUEST)    
